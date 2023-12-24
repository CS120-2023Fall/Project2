#pragma once
#include<deque>
#include<assert.h>
#include<JuceHeader.h>
#include<chrono>
#include<cstdlib>
#include "receiver_transfer.h"
#include "macros.h"

// milisecond
#define ACK_TIME_OUT_THRESHOLD 1200
#define RECEND_THRESHOLD 4

class MAC_Layer {
public:
    MAC_Layer() {}
    MAC_Layer(juce::Label *labels[], int num_labels) {
        if (num_labels > 5) {
            assert(0);
        }
        for (int i = 0; i < num_labels; ++i) {
            mes[i] = labels[i];
        }
    };

    ~MAC_Layer() {
    }
    // update MAC states
    void refresh_MAC(const float *inBuffer, float *outBuffer, int num_samples);
    // prepare for next packet
    void Start() {
        macState = MAC_States_Set::Idle;
        receiver.Initialize();
        transmitter.Initialize();
        resend = 0;
        ackTimeOut_valid = false;
        wait = false;
        backoff_exp = 0;
        startTransmitting = START_TRANS_FIRST;
    }
    
    //void reset_receiving_info();
    void STOP() {
        receiver.Write_symbols();       
    }

public:
    enum class MAC_States_Set {
        Idle,
        CarrierSense,
        RxFrame,
        TxFrame,
        TxACK,
        ACKTimeout,
        LinkError,
        debug_error
    };

    MAC_States_Set macState{MAC_States_Set::Idle};
    bool TxPending{ false };
    std::deque<int> received_data;    
    bool wait = false; // wait for ack
    int start_for_wait_sample=0;
    bool startTransmitting{ START_TRANS_FIRST };

private:
    int mac_address{ MY_MAC_ADDRESS };
    juce::Label *mes[5]{ nullptr }; // array of pointers to send message    
    int resend{ 0 }; // the number of resending times
    // ack time out detect
    // std::chrono::steady_clock::now()
    std::chrono::time_point<std::chrono::steady_clock> beforeTime_ack;
    bool ackTimeOut_valid{ false };    
    int backoff_exp{ 1 }; // exponent of the backoff. 2^m - 1, millisecond
    std::chrono::time_point < std::chrono::steady_clock> beforeTime_backoff{ std::chrono::steady_clock::now() };
public:
    Receiver receiver;
    Transfer transmitter;
};
void KeepSilence(const float* inBuffer, float* outBuffer, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        outBuffer[i] = 0;
    }
}
void MAC_Layer::refresh_MAC(const float *inBuffer, float *outBuffer, int num_samples) {
    /// Idle
    if (macState == MAC_States_Set::Idle) {
        auto currentTime = std::chrono::steady_clock::now();
        // 1. ack time out
        // ///////////////////////
        // pass ackTimeout state, exit directly
        ///////////////////////////////
        if (ackTimeOut_valid) {
            // milisecond
            double duration_millsecond = std::chrono::duration<double, std::milli>(currentTime - beforeTime_ack).count();
            if (duration_millsecond > ACK_TIME_OUT_THRESHOLD) {
                macState = MAC_States_Set::ACKTimeout;//resend the package
                ackTimeOut_valid = false;
                /////////////////////////////// watch out here!!!!!!!!!! ///////////////
                //macState = MAC_States_Set::LinkError;
                /////////////////////////
                return;
            }
        }
        // 2. detect preamble, invoke detect_frame()
        bool tmp = receiver.detect_frame(inBuffer, outBuffer, num_samples);
        if (tmp) {
            mes[3]->setText("preamble detecked " + std::to_string(receiver.received_packet) + ", " + std::to_string(transmitter.transmitted_packet), 
                juce::NotificationType::dontSendNotification);
            macState = MAC_States_Set::RxFrame;
            std::cout << "detect_frame" << std::endl;
            //std::cout <<"receiver_buffer:" <<receiver.receive_buffer.size() << std::endl;
            std::cout << "sync_buffer:" << receiver.sync_buffer.size() << std::endl;
            std::cout << "decode_buffer:" << receiver.decode_buffer.size() << std::endl;
            std::cout << "symbol_code:" << receiver.symbol_code.size() << std::endl;
            // It must return due to the implementation of detect_frame().
            return;
        }
        // 3. send data
        double duration_milisecond = std::chrono::duration<double, std::milli>(currentTime - beforeTime_backoff).count();
        // +, - are prior to <<
        double backoff = (1 << backoff_exp) - 1;
        if (!CSMA_ONLY_RECEIVE && TxPending && (backoff == 0 || duration_milisecond >= backoff)) {
            backoff_exp = 0;
            macState = MAC_States_Set::CarrierSense;
        }
    }
    /// RxFrame
    if (macState == MAC_States_Set::RxFrame) {
        Rx_Frame_Received_Type tmp = receiver.decode_one_packet(inBuffer, outBuffer, num_samples);

        std::cout << "received packet type: " << (int)tmp << std::endl;
        switch (tmp) {
            case Rx_Frame_Received_Type::error:
                macState = MAC_States_Set::Idle;
                return;
            case Rx_Frame_Received_Type::still_receiving:
                return;
            case Rx_Frame_Received_Type::valid_ack:
                ackTimeOut_valid = false;
                transmitter.transmitted_packet += 1;//the next staus transmit the next packet
                macState = MAC_States_Set::Idle;
                mes[2]->setText("Received ack, transmitted packet: " + std::to_string(transmitter.transmitted_packet), 
                    juce::NotificationType::dontSendNotification);
                wait = false;
                backoff_exp = rand() % 5 + 4;
                return;
            case Rx_Frame_Received_Type::valid_data:
                //std::cout << "receiver_buffer:" << receiver.receive_buffer.size() << std::endl;
                std::cout << "sync_buffer:" << receiver.sync_buffer.size() << std::endl;
                std::cout << "decode_buffer:" << receiver.decode_buffer.size() << std::endl;
                std::cout << "symbol_code:" << receiver.symbol_code.size() << std::endl;
                macState = MAC_States_Set::TxACK;
                receiver.received_packet += 1;
                bool feedback = transmitter.Add_one_packet(inBuffer, outBuffer, num_samples, 
                    Tx_frame_status::Tx_ack, receiver.received_packet);
                //backoff_exp = rand() % 5 + 3;
                //beforeTime_backoff = std::chrono::steady_clock::now();
                mes[1]->setText("Packet received: " + std::to_string(receiver.received_packet), juce::dontSendNotification);
                /////////////////////// delete me ////////////////////
                //if (receiver.received_packet * NUM_PACKET_DATA_BITS >= 50000 && 
                //    transmitter.transmitted_packet * NUM_PACKET_DATA_BITS >= 50000 ) {
                //    macState = MAC_States_Set::LinkError;
                //    std::cout << "transmit finish" << std::endl;
                //}
                //////////////////////////////////////////////////////////
                break;
        }
    }
    /// TxACK
    if (macState == MAC_States_Set::TxACK) {
        
        auto currentTime = std::chrono::steady_clock::now();
        double duration_milisecond = std::chrono::duration<double, std::milli>(currentTime - beforeTime_backoff).count();
        // +, - first, then << so use ()
        double backoff = (1 << backoff_exp) - 1;
        if (duration_milisecond <= backoff) {
            return;
        }
        std::cout << "sending ack" << std::endl;
        //if (!receiver.if_channel_quiet(inBuffer, num_samples)) {
        //    return;
        //}
        bool finish = transmitter.Trans(inBuffer, outBuffer, num_samples);
        if (finish) {
            backoff_exp = rand() % 5 + 4;
            macState = MAC_States_Set::Idle;
        }
        // The computer has received a packet. It can start to transmit.
        startTransmitting = true;
        return;
    }
    /// CarrierSense
    if (macState == MAC_States_Set::CarrierSense) {
        if (receiver.if_channel_quiet(inBuffer, num_samples)) {
            macState = MAC_States_Set::TxFrame;
            bool feedback = transmitter.Add_one_packet(inBuffer, outBuffer, num_samples, Tx_frame_status::Tx_data);
        }
        else {
            backoff_exp = rand() % 5 + 4;
            beforeTime_backoff = std::chrono::steady_clock::now();
            macState = MAC_States_Set::Idle;
            return;
        }
    }
    /// TxFrame
    if (macState == MAC_States_Set::TxFrame) {
        bool finish= transmitter.Trans(inBuffer, outBuffer, num_samples);
        // transmition finishes
        if (finish) {
            beforeTime_ack = std::chrono::steady_clock::now();
            ackTimeOut_valid = true;
            macState = MAC_States_Set::Idle;
            wait = true;
        }
        return;
    }
    /// ACKTimeout
    if (macState == MAC_States_Set::ACKTimeout) {
        if (resend > RECEND_THRESHOLD) {
            macState = MAC_States_Set::LinkError;
        }
        // should not reach here
        else {
            ++resend;
            // set backoff after ack timeout
            // [3, 8]
            backoff_exp = rand() % 6 + 3;
            beforeTime_backoff = std::chrono::steady_clock::now();
            macState = MAC_States_Set::Idle;
            wait = false;
            return;
        }
    }
    /// LinkError
    if (macState == MAC_States_Set::LinkError) {
        return;
    }
}
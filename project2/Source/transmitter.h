#ifndef TRANSMITTER_H
#define TRANSMITTER_H
#include <vector>
#include <string>
#include <deque>
#include <iostream>
#include "CRC.h"
#include <cmath>
#include <JuceHeader.h>
#include <iomanip>  // Includes ::std::hex
#include <iostream> // Includes ::std::cout
#include <cstdint>  // Includes ::std::uint32_t
#include "utility.h"

// These macros are aborted.
/////////////////////////////////////////////////////////////
//#define MAXIMUM_TOTAL_SECONDS 5//total time it will run most
//#define PREAMBLE_SIZE 512 // the number of preamble samples
//#define PACKET_DATA_SIZE 500 ////the number of packet, not used by fbz!!!
// I think PACKET_DATA_SIZE is the number of symbols in a packet.
// --fbz
//#define END_FRAME_SIZE 48 //the number of empty zero between packets
//const double pi = juce::MathConstants<double>::pi;
//bool psk_mode = true;
////fsk_mode
//bool test_for_repeat_preamble = false;
////
//bool up_down = true;
//constexpr const double freq_carrier_0 = 3000;
//constexpr const double freq_carrier_1 = 5000;
////psk_mode
//constexpr const double freq_carrier_psk = 10000;
//constexpr const double freq_up_preamble = 1000;
//constexpr const double freq_down_preamble = 1000;
//constexpr const int default_sample_rate = 48000;
//constexpr const int samples_per_bit = default_sample_rate / 1000 / 2;
// TRANSMITTER_CONSTANTS
/////////////////////////////////////////////////////////////////////

// This class is aborted.

//class Transmitter {
//public:
//    Transmitter(const std::string &path, int _sample_rate) :sample_rate(_sample_rate) {
//        bits = Read_bits(path);
//        generate_length();
//        generate_crc();
//    }
//    void generate_crc() {
//        int size;
//        unsigned char *data = unsigned_int_to_unsigned_char_star(from_bits_vector_to_unsigned_int(bits), size);
//        CRC_8 = CRC::CalculateBits(data, size * 8, CRC::CRC_8());
//    }
//
//    void generate_length() {
//        unsigned int size = bits.size();
//        //the maximum bits is 2^15
//        const int maximum_length_bits = 15;
//
//        int l = 1 << (maximum_length_bits - 1);
//        int current_shift = maximum_length_bits - 1;
//
//        while (l) {
//            bool bit = (bool)((l & size) >> current_shift);
//            l = l >> 1;
//            current_shift -= 1;
//            length_bits.push_back(bit);
//        }
//
//    }
//    std::vector<bool> bits;
//    std::vector<double>preamble;
//    std::vector<double> preamble_reverse;
//    std::vector<double>time;
//    std::vector<double> carrier_waves_0;
//    uint8_t CRC_8;
//    std::vector<double> carrier_waves_1;
//    std::vector<bool> length_bits;
//    std::vector<bool> CRC_bits;
//    int sample_rate;
//    std::vector<double> packet_sequences;
//
//};
// a transmitter_with_psk_ask

//#define BITS_PER_SYMBOL 4// the symbol represent how many bit
//float max_amplitude = 1;
//bool check_crc = false;//the crc is 8 bit 
//constexpr const unsigned int CRC_SYMBOL_SIZE = 8 / BITS_PER_SYMBOL;
  //!CRC_CONSTANTS


class Transmitter_with_wire {
public:
    Transmitter_with_wire() = default;
    Transmitter_with_wire(const std::string &path, int _sample_rate) :sample_rate(_sample_rate) {
        bits = Read_bits_from_bin(path);//read the bits from a bin file
        // every 500 bits have a crc calculation
        for (int i = 0; i < bits.size(); i += 500) {
            char tmp = 0;
            for (int j = 0; j < 496; ++j) {
                tmp = (tmp << 1) + (char)bits[i + j];
                if ((j + 1) % 8 == 0) {
                    bytes_for_crc_calculation[j / 8] = tmp;
                    tmp = 0;
                }
            }
            for (int j = 496; j < 500; ++j) {
                tmp = (tmp << 1) + bits[i + j];
            }
            bytes_for_crc_calculation[62] = tmp;
            std::uint32_t crc = CRC::CalculateBits(bytes_for_crc_calculation, sizeof(bytes_for_crc_calculation) - 4, CRC::CRC_32());
            crc_32_t.emplace_back((int)crc);
        }
        generate_preamble();
    }

    void generate_preamble() {
        preamble = std::vector<int>(64, 0);
        // 10101010 * 7 + 10101011
        for (int i = 0; i <= 62; ++i) {
            preamble[i] = (i + 1) % 2;
        }
        preamble[63] = 1;
    }

    // 500 date bits have 32-bit crc.
    void generate_crc() {
        CRC::Table<std::uint32_t, 32> crc_table(CRC::CRC_32());
        //std::uint32_t crc = CRC::Calculate()
    }


    //void generate_crc() {
    //    int size = bits.size();
    //    unsigned int data_bits_per_packet = BITS_PER_SYMBOL * PACKET_DATA_SIZE;
    //    for (int i = 0; i < size; i += data_bits_per_packet) {
    //        std::vector<bool>bits_slice;
    //        uint8_t crc_8;
    //        for (int j = i; j < i + data_bits_per_packet; j++) {
    //            bits_slice.push_back(bits[j]);
    //        }
    //        int data_size;
    //        unsigned char *data = unsigned_int_to_unsigned_char_star(from_bits_vector_to_unsigned_int(bits_slice), data_size);
    //        crc_8 = CRC::CalculateBits(data, data_bits_per_packet, CRC::CRC_8());
    //        std::vector<bool> crc_8_vector = from_uint8_t_to_bits_vector(crc_8);
    //        for (int i = 0; i <= 7; i++) {
    //            CRC_bits.push_back(crc_8_vector[i]);
    //        }
    //    }
    //}

    //void generate_knots() {
    //    unsigned int ask_num = 1 << (seperation_num - 1);
    //    amplitude_step = max_amplitude / (ask_num);
    //}
    //void Write_symbols() {
    //    FILE *file = fopen("symbols.txt", "w");
    //    for (auto &symbol : symbols) {
    //        fprintf(file, "%d ", symbol >> 1);
    //    }
    //    fclose(file);
    //}
    //std::vector<double> generate_one_packet_modulation() {
    //}
    //void generate_crc_one_packet() {
    //}

    //void generate_length() {
    //    unsigned int size = bits.size();
    //    //the maximum bits is 2^15
    //    const int maximum_length_bits = 15;
    //    int l = 1 << (maximum_length_bits - 1);
    //    int current_shift = maximum_length_bits - 1;
    //    while (l) {
    //        bool bit = (bool)((l & size) >> current_shift);
    //        l = l >> 1;
    //        current_shift -= 1;
    //        length_bits.push_back(bit);
    //    }
    //}
    //void generate_carrier() {
    //    double step = 1.0 / static_cast<double>(sample_rate);//the total time step dt per sample
    //    for (int i = 0; i < sample_rate * MAXIMUM_TOTAL_SECONDS; i++) {
    //        time.push_back(i * step);
    //    }
    //    for (int i = 0; i < time.size(); i++) {
    //        //choose the psk modeluation mode
    //        carrier_waves_0.push_back(cos(2 * pi * freq_carrier_psk * time[i] + pi));
    //        carrier_waves_1.push_back(cos(2 * pi * freq_carrier_psk * time[i]));
    //    }
    //}//generate the carrier_waves
    //unsigned int demodulate(double A) {//given the amplitude find the closest to match 
    //    //translate the waves to unsigned int
    //    unsigned int ask_num = 1 << (seperation_num - 1);
    //    unsigned int index = ask_num - 1;
    //    //first find the nearest
    //    float nearest = 0;
    //    for (int i = 0; i < ask_num; i++) {
    //        if (A > i * amplitude_step && A <= (i + 1) * amplitude_step) {
    //            if (A - i * amplitude_step < (i + 1) * amplitude_step - A) {
    //                nearest = i * amplitude_step;
    //                index = i - 1;
    //                if (i == 0) {
    //                    index = 0;
    //                }
    //            }
    //            else {
    //                nearest = (i + 1) * amplitude_step;
    //                index = i;
    //            }
    //            break;
    //        }
    //    }
    //    return index;
    //}
    //void generate_packet_sequences() {// add the head and modulate the data
    //    //the total sequence 
    //    int size = symbols.size();
    //    unsigned int ask_num = 1 << (seperation_num - 1);// if 4 bits then 8 amplitude
    //    for (int i = 0; i < size; i += PACKET_DATA_SIZE) {
    //        int start = i;
    //        int end = i + PACKET_DATA_SIZE;
    //        if (end >= size) {
    //            end = size - 1;
    //        }
    //        for (int j = 0; j < PREAMBLE_SIZE; j++) {
    //            packet_sequences.push_back(preamble[j]);
    //        }
    //        // push the preamble
    //      //then the crc
    //        for (int symbol_index = start; symbol_index <= end; symbol_index++) {
    //            unsigned int symbol = symbols[symbol_index];
    //            auto carrier = (symbol & 1) == 1 ? carrier_waves_1 : carrier_waves_0;//using the last bit to determine psk
    //            float amplitude = max_amplitude * float((symbol >> 1) + 1) / ask_num;
    //            for (int j = 0; j < samples_per_bit; j++) {
    //                packet_sequences.push_back(amplitude * carrier[j]);
    //            }
    //        }
    //        //modulation
    //    }
    //}

    std::vector<bool> bits;
    std::vector<int>preamble;
    int transmitted_packet;
    std::deque<double> transmittion_buffer;
    int sample_rate;

    // crc result, 32 bit width data
    std::vector<int> crc_32_t;
    // crc needs char[]
    char bytes_for_crc_calculation[63]{ 0 };
    //int seperation_num;
    //std::vector<unsigned int > symbols;
    //std::vector<double> preamble_reverse;
    //std::vector<double>time;
    //std::vector<double> carrier_waves_0;
    //std::vector<double> carrier_waves_1;
    //std::vector<bool> length_bits;
    //std::vector<double> sequence;
    //std::vector<double> packet_sequences;
    //float amplitude_step;
    //std::vector<double>knots;
};
Transmitter_with_wire default_trans_wire("INPUT.bin", 48000);

//Transmitter default_trans("INPUT.txt", default_sample_rate);


#endif

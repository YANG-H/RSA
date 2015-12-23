#include <chrono>
#include <iostream>
#include <sstream>

#include "rsa.hpp"

using namespace std::chrono;

high_resolution_clock::time_point last_time_point;
void tic() { last_time_point = high_resolution_clock::now(); }
void toc() {
  auto dur = high_resolution_clock::now() - last_time_point;
  std::cout << "time cost: " << duration_cast<duration<double>>(dur).count()
            << "s" << std::endl;
}

int main(int argc, char **argv) {
  static const int SIZE = 64;

  // setup random engine
  std::default_random_engine rnd;
  rnd.seed(uint32_t(high_resolution_clock::now().time_since_epoch().count()));

  rsa::cipher<SIZE> cipher;
  rsa::decipher<SIZE> decpher;

  std::cout << "bit length: " << SIZE * 32 << std::endl;
  std::cout << " ------- generating keys ------- " << std::endl;
  tic();
  rsa::generate_keys<SIZE>(cipher, decpher);
  std::cout << " ------- keys generated ------- " << std::endl;
  toc();
  std::cout << "==============================================" << std::endl;

  bool firsttime = true;
  while (true) {

    std::string msg;
    std::cout << (firsttime ? " ------- input your message: ------- "
                            : " ------- input your next message: ------- ")
              << std::endl;
    ;
    std::getline(std::cin, msg);
    std::cout << std::endl;

    std::cout << " ------- original message ------- " << std::endl;
    std::cout << msg << std::endl;
    std::cout << std::endl;

    // encode string to big numbers
    auto message_number = cipher.encode(msg);
    std::cout << " ------- message number ------- " << std::endl;
    for (auto i : message_number)
      std::cout << i << std::endl;
    std::cout << std::endl;

    // encrypt
    tic();
    auto encrypted_message_number = cipher.encrypt(message_number);
    std::cout << " ------- encrypted message number ------- " << std::endl;
    for (auto i : encrypted_message_number)
      std::cout << i << std::endl;
    toc();
    std::cout << std::endl;

    // decrypt
    tic();
    auto decrypted_message_number = decpher.decrypt(encrypted_message_number);
    std::cout << " ------- decrypted message number ------- " << std::endl;
    for (auto i : decrypted_message_number)
      std::cout << i << std::endl;
    toc();
    std::cout << std::endl;

    // decode big numbers to string
    std::string recovered_msg = decpher.decode(decrypted_message_number);
    std::cout << " ------- recovered message ------- " << std::endl;
    std::cout << recovered_msg << std::endl << std::endl;

    firsttime = false;
  }

  return 0;
}

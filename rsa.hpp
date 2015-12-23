#pragma once

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "big_num.hpp"

namespace rsa {
template <int len> class cipher;
template <int len> class decipher;

template <int len> class cipher {
  static_assert(sizeof(unsigned) == 4, "incorrect architecture!");

public:
  typedef std::vector<fixed_uint<len>> encoded_sequence_t;

  static encoded_sequence_t encode(const std::string &msg);

  encoded_sequence_t encrypt(const encoded_sequence_t &msg) const;
  inline encoded_sequence_t operator()(const std::string &msg) const {
    return encrypt(encode(msg));
  }

private:
  fixed_uint<len> n, e;
  template <int l> friend void generate_keys(cipher<l> &pbk, decipher<l> &pvk);
};

template <int len>
typename cipher<len>::encoded_sequence_t
cipher<len>::encrypt(const encoded_sequence_t &msg) const {
  encoded_sequence_t m(msg);
  for (auto &i : m)
    i = i.mod_pow(e, n);
  return m;
}

template <int len>
typename cipher<len>::encoded_sequence_t
cipher<len>::encode(const std::string &msg) {
  if (msg.empty())
    return encoded_sequence_t();
  encoded_sequence_t encoded(1 + msg.size() / (3 * len) +
                             1); // encode 3 chars into one int
  encoded[0][0] = uint32_t(msg.size());
  int i = 0, j = 1, k = 0;
  for (char c : msg) {
    encoded[i][j] |= (uint32_t(c) << (8 * k));
    k++;
    if (k == 3) { // encode 3 chars into one int
      k = 0;
      j++;
    }
    if (j == len) {
      j = 0;
      i++;
    }
  }
  return encoded;
}

template <int len> class decipher {
  static_assert(sizeof(unsigned) == 4, "incorrect architecture!");

public:
  typedef std::vector<fixed_uint<len>> encoded_sequence_t;
  static std::string decode(const encoded_sequence_t &encoded);

  encoded_sequence_t decrypt(const encoded_sequence_t &msg) const;
  inline std::string operator()(const encoded_sequence_t &encrypted) const {
    return decode(decrypt(encrypted));
  }

private:
  fixed_uint<len> n, e, d, phi, p, q;
  template <int l> friend void generate_keys(cipher<l> &pbk, decipher<l> &pvk);
};

template <int len>
typename decipher<len>::encoded_sequence_t
decipher<len>::decrypt(const encoded_sequence_t &msg) const {
  encoded_sequence_t m(msg);
  for (auto &i : m)
    i = i.mod_pow(d, n);
  return m;
}

template <int len>
std::string decipher<len>::decode(const encoded_sequence_t &encoded) {
  if (encoded.empty())
    return std::string();
  unsigned msgsize = encoded[0][0];
  std::string msg(msgsize, '\0');
  int i = 0, j = 1, k = 0;
  for (char &c : msg) {
    c = char((encoded[i][j] & (0xffu << (8 * k))) >> (8 * k));
    k++;
    if (k == 3) { // encode 3 chars into one int
      k = 0;
      j++;
    }
    if (j == len) {
      j = 0;
      i++;
    }
  }
  return msg;
}

template <int len> void generate_keys(cipher<len> &pbk, decipher<len> &pvk) {
  typedef fixed_uint<len / 2> half_num_t;
  typedef fixed_uint<len> num_t;

  std::default_random_engine rnd;
  rnd.seed(uint32_t(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()));

  half_num_t p, q;
  num_t n, phi, e, d;

  while (true) {
    std::cout << " ------- locating prime numbers p and q ... ------- "
              << std::endl;

    p.set_to_probable_prime(rnd);
    q.set_to_probable_prime(rnd);

    std::cout << "p: " << std::endl << p << std::endl;
    std::cout << "q: " << std::endl << q << std::endl;

    half_num_t p_1 = p, q_1 = q;
    --p_1;
    --q_1;

    std::cout << "p-1: " << std::endl << p_1 << std::endl;
    std::cout << "q-1: " << std::endl << q_1 << std::endl;

    n = mult_unlimit(p, q);
    phi = mult_unlimit(p_1, q_1);

    std::cout << "n = p*q: " << std::endl << n << std::endl;
    std::cout << "phi = (p-1)*(q-1): " << std::endl << phi << std::endl;

    std::cout << " ------- locating e and d ... ------- " << std::endl;
    for (int i = 0; i < 100; i++) {
      e.set_to_probable_prime(rnd, 1);
      if (e >= phi)
        continue;

      d = e.mod_inverse(phi);
      std::cout << "e: " << std::endl << e << std::endl;
      std::cout << "d = e^-1 mod phi: " << std::endl << d << std::endl;
      std::cout << "e*d mod phi: " << std::endl
                << e.mod_mult(d, phi) << std::endl;

      pbk.n = n;
      pbk.e = e;

      pvk.n = n;
      pvk.e = e;
      pvk.d = d;
      pvk.phi = phi;
      pvk.p = num_t(p);
      pvk.q = num_t(q);

      return;
    }
  }
}
}

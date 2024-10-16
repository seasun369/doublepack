/**
 * @file channel.h
 *
 * SCL --- Secure Computation Library
 * Copyright (C) 2022 Anders Dalskov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 * USA
 */
#ifndef _SCL_NET_CHANNEL_H
#define _SCL_NET_CHANNEL_H

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "scl/math/mat.h"
#include "scl/math/vec.h"
#include "scl/net/config.h"
#include <NTL/ZZ_pE.h>
#include <sstream>

namespace scl {

/**
 * @brief The maximum size that scl::Channel will read when receiving Vecs.
 *
 * This is used in order to guard against malicious inputs and to avoid lockups
 * in case of programming mistakes. The default size <code>1 << 25</code>
 * corresponds to reading a vector of roughly 500mb in when field elements take
 * up 16 bytes.
 */
#ifndef MAX_VEC_READ_SIZE
#define MAX_VEC_READ_SIZE 1 << 25
#endif

/**
 * @brief The maxium amount of elements that scl::Mat::Read will read.
 */
#ifndef MAX_MAT_READ_SIZE
#define MAX_MAT_READ_SIZE 1 << 25
#endif

namespace {

template <typename T>
using IsCopyable =
    typename std::enable_if_t<std::is_trivially_copyable_v<T>, bool>;

}  // namespace

#define _SCL_CC(x) reinterpret_cast<const unsigned char*>(x)
#define _SCL_C(x) reinterpret_cast<unsigned char*>(x)

/**
 * @brief Abstract channel for communicating between two peers.
 *
 * scl::Channel defines the interface for a channel between two peers, as well
 * as a number of convenience methods for sending and receiving different kinds
 * of data. To implement an actual channel, subclass scl::Channel and implement
 * the four virtual methods.
 *
 * @see scl::InMemoryChannel
 * @see scl::TcpChannel
 */
class Channel {
 public:
  /**
   * @brief Close connection to remote.
   */
  virtual void Close() = 0;

  /**
   * @brief Send data to the remote party.
   * @param src the data to send
   * @param n the number of bytes to send
   */
  virtual void Send(const unsigned char* src, std::size_t n) = 0;

  /**
   * @brief Receive data from the remote party.
   * @param dst where to store the received data
   * @param n how much data to receive
   */
  virtual void Recv(unsigned char* dst, std::size_t n) = 0;

  /**
   * @brief Send a trivially copyable item.
   * @param src the thing to send
   */
  template <typename T, IsCopyable<T> = true>
  void Send(const T& src) {
    Send(_SCL_CC(&src), sizeof(T));
  }

  /**
   * @brief Send a vector of trivially copyable things.
   *
   * <code>src.size()</code> is used to determine how many bytes to read, so \p
   * src must have enough room for the data we expect to receive.
   *
   * @param src an STL vector of things to send
   */
  template <typename T, IsCopyable<T> = true>
  void Send(const std::vector<T>& src) {
    Send(_SCL_CC(src.data()), sizeof(T) * src.size());
  }

  /**
   * @brief Send a vector object.
   *
   * Note that \p T cannot be guaranteed to be trivially copyable so this method
   * needs to make a temporary copy of \p vec in order to serialize it
   * correctly. If \p T is known to be trivially copyable, then it might be
   * faster to call <code>Send(vec.ToStlVector())</code> instead.
   *
   * This method sends the size of the vector first followed by its content.
   *
   * @param vec the Vector
   */
  template <typename T>
  void Send(const Vec<T>& vec) {
    const std::uint32_t vec_size = vec.Size();
    Send(vec_size);
    // have to make a copy here since it's not guaranteed that we can directly
    // write T to the channel.
    auto buf = std::make_unique<unsigned char[]>(vec.ByteSize());
    vec.Write(buf.get());
    Send(_SCL_CC(buf.get()), vec.ByteSize());
  }

  /**
   * @brief Send a matrix.
   *
   * Like with Send(const Vec<T>&) this method makes a copy of \p mat internally
   * in order to serialize the matrix correctly.
   *
   * Also similar to Send(const Vec<T>&) this method first sends the row count,
   * then column count and finally the matrix content.
   *
   * @param mat the matrix to send
   */
  template <typename T>
  void Send(const Mat<T>& mat) {
    const std::uint32_t rows = mat.Rows();
    const std::uint32_t cols = mat.Cols();
    Send(rows);
    Send(cols);
    auto buf = std::make_unique<unsigned char[]>(mat.ByteSize());
    mat.Write(buf.get());
    Send(_SCL_CC(buf.get()), mat.ByteSize());
  }

  /**
   * @brief Receive a trivially copyable item.
   * @param dst where to store the received item
   */
  template <typename T, IsCopyable<T> = true>
  void Recv(T& dst) {
    Recv(_SCL_C(&dst), sizeof(T));
  }

  /**
   * @brief Receive a vector of trivially copyable items.
   *
   * <code>dst.size()</code> determines how many bytes to receive.
   *
   * @param dst where to store the received items
   */
  template <typename T, IsCopyable<T> = true>
  void Recv(std::vector<T>& dst) {
    Recv(_SCL_C(dst.data()), sizeof(T) * dst.size());
  }

  /**
   * @brief Receive a vector.
   * @param vec where to store the received vector
   * @throws std::logic_error in case the received vector size exceeds
   * MAX_VEC_READ_SIZE
   */
  template <typename T>
  void Recv(Vec<T>& vec) {
    auto vec_size = RecvSize();
    if (vec_size > MAX_VEC_READ_SIZE)
      throw std::logic_error("received vector exceeds size limit");
    auto n = vec_size * T::ByteSize();
    auto buf = std::make_unique<unsigned char[]>(n);
    Recv(_SCL_C(buf.get()), n);
    vec = Vec<T>::Read(vec_size, _SCL_CC(buf.get()));
  }

  /**
   * @brief Receive a Matrix.
   * @param mat where to store the received matrix
   * @throws std::logic_error in case the row count times column count exceeds
   * MAX_MAT_READ_SIZE
   */
  template <typename T>
  void Recv(Mat<T>& mat) {
    auto rows = RecvSize();
    auto cols = RecvSize();
    if (rows * cols > MAX_MAT_READ_SIZE)
      throw std::logic_error("received matrix exceeds size limit");
    auto n = rows * cols * T::ByteSize();
    auto buf = std::make_unique<unsigned char[]>(n);
    Recv(_SCL_C(buf.get()), n);
    mat = Mat<T>::Read(rows, cols, _SCL_CC(buf.get()));
  }

  // 序列化 ZZ_pE
  inline std::vector<unsigned char> serializeZZ_pE(const NTL::ZZ_pE& value) {
    std::ostringstream oss;
    oss << value;
    std::string str = oss.str();
    return std::vector<unsigned char>(str.begin(), str.end());
  }

  // 反序列化 ZZ_pE
  inline NTL::ZZ_pE deserializeZZ_pE(const std::vector<unsigned char>& buffer) {
    std::string str(buffer.begin(), buffer.end());
    std::istringstream iss(str);
    NTL::ZZ_pE value;
    iss >> value;
    return value;
  }

  // 为 NTL::ZZ_pE 类型添加专门的 Send 和 Recv 函数
  void Send(const NTL::ZZ_pE& src) {
    std::vector<unsigned char> buffer = serializeZZ_pE(src);
    Send(static_cast<std::uint32_t>(buffer.size()));  // 先发送大小
    Send(buffer.data(), buffer.size());  // 再发送数据
  }

  void Recv(NTL::ZZ_pE& dst) {
    std::uint32_t size;
    Recv(size);  // 先接收大小
    std::vector<unsigned char> buffer(size);
    Recv(buffer.data(), size);  // 再接收数据
    dst = deserializeZZ_pE(buffer);
  }

  // 序列化 ZZ_p
  inline std::vector<unsigned char> serializeZZ_p(const NTL::ZZ_p& value) {
    std::ostringstream oss;
    oss << value;
    std::string str = oss.str();
    return std::vector<unsigned char>(str.begin(), str.end());
  }

  // 反序列化 ZZ_p
  inline NTL::ZZ_p deserializeZZ_p(const std::vector<unsigned char>& buffer) {
    std::string str(buffer.begin(), buffer.end());
    std::istringstream iss(str);
    NTL::ZZ_p value;
    iss >> value;
    return value;
  }

  // 为 NTL::ZZ_p 类型添加专门的 Send 函数
  void Send(const NTL::ZZ_p& src) {
    std::vector<unsigned char> buffer = serializeZZ_p(src);
    Send(static_cast<std::uint32_t>(buffer.size()));  // 先发送大小
    Send(buffer.data(), buffer.size());  // 再发送数据
  }

  // 为 NTL::ZZ_p 类型添加专门的 Recv 函数
  void Recv(NTL::ZZ_p& dst) {
    std::uint32_t size;
    Recv(size);  // 先接收大小
    std::vector<unsigned char> buffer(size);
    Recv(buffer.data(), size);  // 再接收数据
    dst = deserializeZZ_p(buffer);
  }

 private:
  std::uint32_t RecvSize() {
    std::uint32_t size;
    Recv(size);
    return size;
  }
};

#undef _SCL_C
#undef _SCL_CC

}  // namespace scl

#endif /* _SCL_NET_CHANNEL_H */

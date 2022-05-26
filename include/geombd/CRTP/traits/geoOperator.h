// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Alvaro Paz

#pragma once

#ifndef GEOMETRIC_OPERATOR_H
#define GEOMETRIC_OPERATOR_H

namespace Eigen {


//template<typename Derived>
//class adSzBase : public ReturnByValue<Derived>
//{
//  private:
//    typedef typename internal::traits<Derived> Traits;
//    typedef typename Traits::Scalar Scalar;

//  protected:
//    typedef typename Traits::Lhs Lhs;
//    typedef typename Traits::Rhs Rhs;

//  public:
//    /*! \brief Constructor. */
//    adSzBase(const Lhs& A, const Rhs& B)
//      : m_A(A), m_B(B)
//    {}

//    inline Index rows() const { return m_A.rows(); }
//    inline Index cols() const { return m_A.cols(); }

//    /*!
//     * This overrides ReturnByValue::coeff because this function is
//     * efficient enough.
//     */
//    Scalar coeff(Index row, Index col) const
//    {
//      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
//             m_B.coeff(row % m_B.rows(), col % m_B.cols());
//    }

//    /*!
//     * This overrides ReturnByValue::coeff because this function is
//     * efficient enough.
//     */
//    Scalar coeff(Index i) const
//    {
//      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
//      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
//    }

//  protected:
//    typename Lhs::Nested m_A;
//    typename Rhs::Nested m_B;
//};

////!-------------------------------------------------------------------------

//template<typename Lhs, typename Rhs>
//class adSz : public adSzBase< adSz<Lhs,Rhs> >
//{
//  private:
//    typedef adSzBase< adSz > Base;
//    using Base::m_A;
//    using Base::m_B;

//  public:
//    /*! \brief Constructor. */
//    adSz(const Lhs& A, const Rhs& B)
//      : Base(A, B)
//    {}

//    /*! \brief Evaluate the Vanish tensor product. */
//    template<typename Dest> void evalTo(Dest& dst) const;
//};

////!-------------------------------------------------------------------------

//template<typename Lhs, typename Rhs>
//template<typename Dest>
//void adSz<Lhs,Rhs>::evalTo(Dest& dst) const
//{
////  const int BlockRows = Rhs::RowsAtCompileTime,
////            BlockCols = Rhs::ColsAtCompileTime;
////  const Index Br = m_B.rows(),
////              Bc = m_B.cols();

//  for (Index j=0; j < m_A.cols(); ++j) dst(0,j) = m_A.coeff(1,j);

//  for (Index j=0; j < m_A.cols(); ++j) dst(1,j) = -m_A.coeff(0,j);

//  for (Index j=0; j < m_A.cols(); ++j) dst(3,j) = m_A.coeff(4,j);

//  for (Index j=0; j < m_A.cols(); ++j) dst(4,j) = -m_A.coeff(3,j);

//}




////    D_q_c_.template row(0) =  D_c_aux.template row(1);  // -ad(Sz)* effect
////    D_q_c_.template row(1) = -D_c_aux.template row(0);
////    D_q_c_.template row(3) =  D_c_aux.template row(4);
////    D_q_c_.template row(4) = -D_c_aux.template row(3);





////!-------------------------------------------------------------------------

//namespace internal {

//template<typename _Lhs, typename _Rhs>
//struct traits< adSz<_Lhs,_Rhs> >
//{
//  typedef typename remove_all<_Lhs>::type Lhs;
//  typedef typename remove_all<_Rhs>::type Rhs;
//  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
//  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

//  enum {
//    Rows = 6,//size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
//    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
//    MaxRows = 6,//size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
//    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
//  };

//  typedef Matrix<Scalar,Rows,Cols> ReturnType;
//};

//} // end namespace internal

////!-------------------------------------------------------------------------

//template<typename A, typename B>
//adSz<A,B> op_adSz(const MatrixBase<A>& a, const MatrixBase<B>& b)
//{
//  return adSz<A, B>(a.derived(), b.derived());
//}




////!-------------------------------------------------------------------------
////!------------------------------NEW CLASS----------------------------------
////!-------------------------------------------------------------------------


//template<typename Derived>
//class Add2ColsBase : public ReturnByValue<Derived>
//{
//  private:
//    typedef typename internal::traits<Derived> Traits;
//    typedef typename Traits::Scalar Scalar;

//  protected:
//    typedef typename Traits::Lhs Lhs;
//    typedef typename Traits::Rhs Rhs;

//  public:
//    /*! \brief Constructor. */
//    Add2ColsBase(const Lhs& A, const Rhs& B)
//      : m_A(A), m_B(B)
//    {}

//    inline Index rows() const { return 6; }
//    inline Index cols() const { return m_A.cols(); }

//    /*!
//     * This overrides ReturnByValue::coeff because this function is
//     * efficient enough.
//     */
//    Scalar coeff(Index row, Index col) const
//    {
//      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
//             m_B.coeff(row % m_B.rows(), col % m_B.cols());
//    }

//    /*!
//     * This overrides ReturnByValue::coeff because this function is
//     * efficient enough.
//     */
//    Scalar coeff(Index i) const
//    {
//      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
//      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
//    }

//  protected:
//    typename Lhs::Nested m_A;
//    typename Rhs::Nested m_B;
//};

////!-------------------------------------------------------------------------

//template<typename Lhs, typename Rhs>
//class Add2Cols : public Add2ColsBase<Add2Cols<Lhs,Rhs> >
//{
//  private:
//    typedef Add2ColsBase<Add2Cols> Base;
//    using Base::m_A;
//    using Base::m_B;

//  public:
//    /*! \brief Constructor. */
//    Add2Cols(const Lhs& A, const Rhs& B)
//      : Base(A, B)
//    {}

//    /*! \brief Evaluate the Kronecker tensor product. */
//    template<typename Dest> void evalTo(Dest& dst) const;
//};

////!-------------------------------------------------------------------------

//template<typename Lhs, typename Rhs>
//template<typename Dest>
//void Add2Cols<Lhs,Rhs>::evalTo(Dest& dst) const
//{
//  const int BlockRows = Rhs::RowsAtCompileTime,
//            BlockCols = Rhs::ColsAtCompileTime;
//  const Index Ar = m_A.cols(),
//              Br = m_B.cols();

////  dst.setZero();

//  for (Index j=0; j < Br; j++) {
//      Block<Dest,6,1>(dst,0,m_B.coeff(j),6,1) = m_A.col(j);
//    }

////  Block<Dest,6,1>(dst,0,j,6,1) = m_A.col(z) + m_A.col(y);
////    for (Index j = Vanisher_i; j < m_A.cols(); j += Vanisher_n)
////      Block<Dest,BlockRows,BlockCols>(dst,i*Br,j*Bc,Br,Bc) = m_A.coeff(i,j) * m_B;

//}

////!-------------------------------------------------------------------------

//namespace internal {

//template<typename _Lhs, typename _Rhs>
//struct traits<Add2Cols<_Lhs,_Rhs> >
//{
//  typedef typename remove_all<_Lhs>::type Lhs;
//  typedef typename remove_all<_Rhs>::type Rhs;
//  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
//  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

//  enum {
//    Rows = 6,//size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
//    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
//    MaxRows = 6,//size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
//    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
//  };

//  typedef Matrix<Scalar,Rows,Cols> ReturnType;
//};

//} // end namespace internal

////!-------------------------------------------------------------------------

//template<typename A, typename B>
//Add2Cols<A,B> add2Cols(const MatrixBase<A>& a, const MatrixBase<B>& b)
//{
//  return Add2Cols<A, B>(a.derived(), b.derived());
//}



////!-------------------------------------------------------------------------
////!-----------------------New Class for Eigen::All--------------------------
////!-------------------------------------------------------------------------


template<typename Derived>
class ColAssignBase : public ReturnByValue<Derived>
{
  private:
    typedef typename internal::traits<Derived> Traits;
    typedef typename Traits::Scalar Scalar;

  protected:
    typedef typename Traits::Lhs Lhs;
    typedef typename Traits::Rhs Rhs;

  public:
    /*! \brief Constructor. */
    ColAssignBase(const Lhs& A, const Rhs& B)
      : m_A(A), m_B(B)
    {}

    inline Index rows() const { return m_A.rows() * m_B.rows(); }
    inline Index cols() const { return m_A.cols() * m_B.cols(); }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index row, Index col) const
    {
      return m_A.coeff(row / m_B.rows(), col / m_B.cols()) *
             m_B.coeff(row % m_B.rows(), col % m_B.cols());
    }

    /*!
     * This overrides ReturnByValue::coeff because this function is
     * efficient enough.
     */
    Scalar coeff(Index i) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
      return m_A.coeff(i / m_A.size()) * m_B.coeff(i % m_A.size());
    }

  protected:
    typename Lhs::Nested m_A;
    typename Rhs::Nested m_B;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
class ColAssign : public ColAssignBase<ColAssign<Lhs,Rhs> >
{
  private:
    typedef ColAssignBase<ColAssign> Base;
    using Base::m_A;
    using Base::m_B;

  public:
    /*! \brief Constructor. */
    ColAssign(const Lhs& A, const Rhs& B)
      : Base(A, B)
    {}

    /*! \brief Evaluate the Kronecker tensor product. */
    template<typename Dest> void evalTo(Dest& dst) const;
};

//!-------------------------------------------------------------------------

template<typename Lhs, typename Rhs>
template<typename Dest>
void ColAssign<Lhs,Rhs>::evalTo(Dest& dst) const
{
  const int BlockRows = Rhs::RowsAtCompileTime,
            BlockCols = Rhs::ColsAtCompileTime;
  const Index Ac = m_A.cols(),
              Bc = m_B.cols();

//  for (int j=0; j < Bc; j++) {
//      dst.col( m_B(j) ) = m_A.col(j);
//    }

  for (int j=0; j < Bc; j++) {
      Block<Dest,6,1>(dst, 0, m_B(j), 6, 1) = m_A.col(j);
    }

}

//!-------------------------------------------------------------------------

namespace internal {

template<typename _Lhs, typename _Rhs>
struct traits<ColAssign<_Lhs,_Rhs> >
{
  typedef typename remove_all<_Lhs>::type Lhs;
  typedef typename remove_all<_Rhs>::type Rhs;
//  typedef typename ScalarBinaryOpTraits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType Scalar;
  typedef double Scalar;
  typedef typename promote_index_type<typename Lhs::StorageIndex, typename Rhs::StorageIndex>::type StorageIndex;

  enum {
    Rows = size_at_compile_time<traits<Lhs>::RowsAtCompileTime, traits<Rhs>::RowsAtCompileTime>::ret,
    Cols = size_at_compile_time<traits<Lhs>::ColsAtCompileTime, traits<Rhs>::ColsAtCompileTime>::ret,
    MaxRows = size_at_compile_time<traits<Lhs>::MaxRowsAtCompileTime, traits<Rhs>::MaxRowsAtCompileTime>::ret,
    MaxCols = size_at_compile_time<traits<Lhs>::MaxColsAtCompileTime, traits<Rhs>::MaxColsAtCompileTime>::ret
  };

  typedef Matrix<Scalar,Rows,Cols> ReturnType;
};

} // end namespace internal

//!-------------------------------------------------------------------------

template<typename A, typename B>
ColAssign<A,B> colAssign(const MatrixBase<A>& a, const MatrixBase<B>& b)
{
  return ColAssign<A, B>(a.derived(), b.derived());
}



} // end namespace Eigen

#endif // GEOMETRIC_OPERATOR_H

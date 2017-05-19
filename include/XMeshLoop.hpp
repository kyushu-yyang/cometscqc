/**
 *  @file   XMeshLoop.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   19 Aug 2016
 **/

#ifndef XMeshLoop_HH
#define XMeshLoop_HH

#include "XCoilBase.hpp"

/// class to handle the 3 dimensional mesh loop
//
class XMeshLoop
{
  public:
    /// @brief  constructor
    /// @detail usage of using the explicit constructor with over one arguments
    ///         \code
    ///          XMeshLoop loop = {mz, mp, mr};
    ///         \endcode
    XMeshLoop();
    explicit XMeshLoop(const int mz, const int mp, const int mr);

    /// @brief deconstructor
    virtual ~XMeshLoop();

    /// @brief initialize the mesh loop
    void SetMesh(const int mz, const int mp, const int mr);

    /// @brief return the mesh
    size_t GetMesh(const Coil dim) const;

    /// @brief check the id is over range or not
    /// @param id the id number of the data
    bool IsOverRange(const int id) const;

    /// @brief  return the id number of 3D data
    /// @param  i id number along the z-axis
    /// @param  j id number along the \f$\phi\f$-axis
    /// @param  k id number along the r-axis
    /// @detail the sequence of loop: \f$ z \rightarrow \phi \rightarrow r \f$.
    //          therefore, id can be calculated as follows: 
    ///         \f$ id = k\cdot (Msh_{z}+2) \cdot (Msh_{\phi}+2) + j\cdot (Msh_{z}+2) + i \f$
    const int Id(const int i, const int j, const int k);

    /// @brief  create a loop for r/phi/z direction
    /// @detail user must implement the function inside the loop
    ///         otherwise, it will loop for nothing
    void StartLoop();

    /// @brief  void function inside the loop
    /// @detail it must be modified by user
    virtual void InLoopR  (int i, int j, int k) {}
    virtual void InLoopPhi(int i, int j, int k) {}
    virtual void InLoopZ  (int i, int j, int k) {}


  protected:
    int fMshZ;
    int fMshP;
    int fMshR;
};

#endif

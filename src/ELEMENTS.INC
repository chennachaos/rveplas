C***********************************************************************
C----------------------------------------------------------------------*
C                                                                      *
C                                                                      *
C                  H Y P L A S   ELEMENTS DATABASE                     *
C                                                                      *
C                                                                      *
C----------------------------------------------------------------------*
C***********************************************************************
C
C
C***********************************************************************
C
C 1. Dimensioning parameters associated with the available elements
C
C***********************************************************************
C
      PARAMETER( MNODE =20          )
      PARAMETER( MDOFN =3           )
      PARAMETER( MTOTG =27          )
      PARAMETER( MEVAB =MNODE*MDOFN )
      PARAMETER( MREPRP=200         )
      PARAMETER( MIEPRP=100         )
C
C
C   * MNODE  = Maximun number of nodes per element.
C
C   * MDOFN  = Maximun number of degrees of freedom per node.
C
C   * MTOTG  = Maximun number of Gauss points allowed in any element.
C
C   * MEVAB  = Maximun number of element variables (degrees of freedom).
C
C   * MREPRP = Maximun number of real element properties that can be
C              stored in array RELPRP for each element type.
C
C   * MIEPRP = Maximun number of integer element properties that can be
C              stored in array IELPRP for each element type.
C
C
C
C***********************************************************************
C
C 2. Element classes enumeration
C
C***********************************************************************
C
      INTEGER  FBAR  ,STDARD
      PARAMETER( STDARD=100 )
      PARAMETER( FBAR  =101 )
C
C
C   * STDARD = Standard displacement based isoparametric elements.
C
C   * FBAR   = F-bar type elements.
C
C
C
C***********************************************************************
C
C 3. Element types enumeration
C
C***********************************************************************
C
      INTEGER  TRI3, TRI6, QUAD4 ,QUAD8 ,QUA4FB ,HEXA8 ,HEX8FB,TETA7
      PARAMETER( TRI3  =1 )
      PARAMETER( TRI6  =7 )
      PARAMETER( QUAD4 =2 )
      PARAMETER( QUAD8 =3 )
      PARAMETER( QUA4FB=4 )
      PARAMETER( HEXA8 =5 )
      PARAMETER( HEX8FB=6 )
      PARAMETER( TETA7 =10)
C
C   * TRI3   = Standard linear 3-noded triangle.
C
C   * TRI6   = Standard quadratic 6-noded triangle.
C
C   * QUAD4  = Standard bi-linear 4-noded quadrilateral.
C
C   * QUAD8  = Standard quadratic 8-noded quadrilateral.
C
C   * QUA4FB = F-Bar 4-noded quadrilateral.
C
C   * HEXA8  = Standard tri-linear 8-noded hexahedron.
C
C   * HEX8FB = F-Bar 8-noded hexahedron.
C
C   * TETA7 = Standard 10-noded Tetrahedron.
C
C***********************************************************************

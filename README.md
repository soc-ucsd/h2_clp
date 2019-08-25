# H2 controller synthesis using closed-loop parameterizations

## Related paper
Y. Zheng, L. Fuerier, M. Kamgarpour, N Li. On the closed-loop parameterization of stabilizing controllers, in final preparation.

## Problem description
 The plant dynamics are:
                           x = Ax_t + Bu_t
                           y = Cx_t + w_t
 The orginal problem is as follows
               min_{K} max_{G} ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
                  subject to     K internally stabilizes G
                                 K \in S
 where Y, U denote the closed-loop transfer matrices from w to y and u, respectively.

 We solve the problem using closed-loop parameterizations, one of them is as follows

              min_{Y,U,W,Z}  ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
               subjec to      [I -G][Y W]
                                    [U Z]  = [I 0]   (1)
                                 [Y W][-G] = [I]
                                 [U Z][I]    [0]     (2)
                             Y,U,W,Z \in FIR(N)     (3)
                              Y \in R,  U \in T      (4)


 Rely on YALMIP to reformulate the above problem into an SDP, then call  Mosek/SeDuMi to get a soluton

## Syntax

[K,H2,info] = clph2(A,B,C,Q,R,userOpts)

 Input variables
      (A,B,C):    system dynamics in discrete time
      Q:    performance weights on output y
      R:    performance weights on input u

 userOpts is a structure and contains the following options
      N:      Oder of FIR approximation    (default:8)
      solver: sedumi, sdpt3, csdp or mosek (default)
      spa:    Distributed control Yes/No   (default: 0)
      S:      Sparsity pattern for the controller  (default: [])

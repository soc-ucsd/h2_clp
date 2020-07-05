# H2 controller synthesis using closed-loop parameterizations

## Problem description
The plant dynamics are:
 
                           x = Ax_t + B(u_t + v_t)
                           y = Cx_t + w_t
                           
The orginal problem is as follows
 
                      min_{K}    ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
                                      + ||Q^{1/2}W||^2 + ||R^{1/2}Z||^2
                  subject to     K internally stabilizes G
                                 K \in S
                                 
where Y, U denote the closed-loop transfer matrices from w to y and u, and W, Z, denote the closed-loop transfer matrices from v to y and u, S is a binary matrix that encodes the controller structure.

 We solve the problem using closed-loop parameterizations, one of them is as follows

              min_{Y,U,W,Z}  ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
                                + ||Q^{1/2}W||^2 + ||R^{1/2}Z||^2
               subjec to      [I -G][Y W]
                                    [U Z]  = [I 0]   (1)
                                 [Y W][-G] = [I]
                                 [U Z][I]    [0]     (2)
                             Y,U,W,Z \in FIR(N)      (3)
                              Y \in R,  U \in T      (4)

where (1)-(3) encode the internal stability constraint, and (4) encodes the sparsity constraint S using the notion of Sparsity invariance.

Rely on YALMIP to reformulate the above problem into an SDP, then call  Mosek/SeDuMi to get a soluton

## Syntax

       >>  [K,H2,info] = clph2(A,B,C,Q,R,userOpts)

1. Input variables
   - (A,B,C):    system dynamics in discrete time
   - Q:    performance weights on output y
   - R:    performance weights on input u
2. userOpts is a structure and contains the following options
   - N:      Order of FIR approximation    (default:8)
   - solver: sedumi, sdpt3, csdp or mosek (default)
   - spa:    Distributed control Yes/No   (default: 0)
   - S:      Sparsity pattern for the controller  (default: [])
      
## Related paper
```
@misc{zheng2019systemlevel,
    title={System-level, Input-output and New Parameterizations of Stabilizing Controllers, and Their Numerical Computation},
    author={Yang Zheng and Luca Furieri and Maryam Kamgarpour and Na Li},
    year={2019},
    eprint={1909.12346},
    archivePrefix={arXiv},
    primaryClass={math.OC}
}
```

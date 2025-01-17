"""
two-site contraction

```
                                            a ────┬──── c          
a ────┬──c ──┬──── f                        │     b     │  
│     b      e     │                        ├─ e ─┼─ f ─┤  
├─ g ─┼─  h ─┼─ i ─┤                        g     h     i 
│     k      n     │                        ├─ j ─┼─ k ─┤ 
j ────┴──l ──┴──── o                        │     m     │ 
                                            l ────┴──── n 
```
"""
oc_H_leg3 = ein"(((agj,abc),gkhb),jkl),(((fio,cef),hnie),lno) -> "
oc_V_leg3 = ein"(((abc,aeg),ehfb),cfi),(gjl,(jmkh,(ikn,lmn))) -> "

"""
two-site contraction

```
                                            a ────┬──── d          
a ────┬──d ──┬──── g                        │     bc    │  
│     bc     ef    │                        ├─ ef─┼─gh ─┤  u
├─ hi─┼─ jk ─┼─lm ─┤                        i     jk    l 
│     op     rs    │                        ├─ mn─┼─op ─┤  v
n ────┴──q ──┴──── t                        │     rs    │ 
      u      v                              q ────┴──── t 
```
"""
oc_H_leg4 = ein"((((ahin,abcd),hojbu),ipkcu),nopq),((((glmt,defg),jrlev),ksmfv),qrst) -> "
oc_V_leg4 = ein"((((abcd,aefi),ejgbu),fkhcu),dghl),(imnq,(mrojv,(nspkv,(lopt,qrst)))) -> "

function contract_n2_H(FLo::leg3, ACu, A1, ACd, FRo, ARu, A2, ARd)
    D1,D2,D3,D4,_ = size(A1)
    M1 = reshape(ein"abcde,fghme->afbgchdm"(A1, conj(A1)), D1^2,D2^2,D3^2,D4^2)
    D1,D2,D3,D4,_ = size(A2)
    M2 = reshape(ein"abcde,fghme->afbgchdm"(A2, conj(A2)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc_H_leg3(FLo, ACu, M1, ACd, FRo, ARu, M2, ARd))
end

function contract_o2_H(FLo::leg3, ACu, A1, ACd, FRo, ARu, A2, ARd, O1, O2)
    D1,D2,D3,D4,_ = size(A1)
    M1 = reshape(ein"(abcde,en),fghmn->afbgchdm"(A1, O1, conj(A1)), D1^2,D2^2,D3^2,D4^2)
    D1,D2,D3,D4,_ = size(A2)
    M2 = reshape(ein"(abcde,en),fghmn->afbgchdm"(A2, O2, conj(A2)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc_H_leg3(FLo, ACu, M1, ACd, FRo, ARu, M2, ARd))
end

function contract_n2_V(ACu::leg3, FLu, A1, FRu, FLo, A2, FRo, ACd)
    D1,D2,D3,D4,_ = size(A1)
    M1 = reshape(ein"abcde,fghme->afbgchdm"(A1, conj(A1)), D1^2,D2^2,D3^2,D4^2)
    D1,D2,D3,D4,_ = size(A2)
    M2 = reshape(ein"abcde,fghme->afbgchdm"(A2, conj(A2)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc_V_leg3(ACu, FLu, M1, FRu, FLo, M2, FRo, ACd))
end

function contract_o2_V(ACu::leg3, FLu, A1, FRu, FLo, A2, FRo, ACd, O1, O2)
    D1,D2,D3,D4,_ = size(A1)
    M1 = reshape(ein"(abcde,en),fghmn->afbgchdm"(A1, O1, conj(A1)), D1^2,D2^2,D3^2,D4^2)
    D1,D2,D3,D4,_ = size(A2)
    M2 = reshape(ein"(abcde,en),fghmn->afbgchdm"(A2, O2, conj(A2)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc_V_leg3(ACu, FLu, M1, FRu, FLo, M2, FRo, ACd))
end

function contract_n2_H(FLo::leg4, ACu, A1, ACd, FRo, ARu, A2, ARd)
    return sum(oc_H_leg4(FLo, ACu, A1, conj(A1), ACd, FRo, ARu, A2, conj(A2), ARd))
end

function contract_o2_H(FLo::leg4, ACu, A1, ACd, FRo, ARu, A2, ARd, O1, O2)
    return sum(oc_H_leg4(FLo, ACu, ein"abcde,ef->abcdf"(A1, O1), conj(A1), ACd, FRo, ARu, ein"abcde,ef->abcdf"(A2, O2), conj(A2), ARd))
end

function contract_n2_V(ACu::leg4, FLu, A1, FRu, FLo, A2, FRo, ACd)
    return sum(oc_V_leg4(ACu, FLu, A1, conj(A1), FRu, FLo, A2, conj(A2), FRo, ACd))
end

function contract_o2_V(ACu::leg4, FLu, A1, FRu, FLo, A2, FRo, ACd, O1, O2)
    return sum(oc_V_leg4(ACu, FLu, ein"abcde,ef->abcdf"(A1, O1), conj(A1), FRu, FLo, ein"abcde,ef->abcdf"(A2, O2), conj(A2), FRo, ACd))
end

"""
one-site contraction

```                         
a ────┬──── c                      
│     b     │                      
├─ e ─┼─ f ─┤                      
│     h     │                      
g ────┴──── i                      
                            
```
"""
oc1_leg3 = ein"(((aeg,abc),ehfb),ghi),cfi -> "

"""
one-site contraction

```                         
a ────┬──── d                      
│     bc    │                      
├─ ef─┼─gh ─┤                      
│     jk    │                      
i ────┴──── l                      
                            
```
"""
oc1_leg4 = ein"((((aefi,abcd),ejgbu),fkhcu),ijkl),dghl -> "

function contract_n1(FLo::leg3, ACu, A, ACd, FRo)
    D1,D2,D3,D4,_ = size(A)
    M = reshape(ein"abcde,fghme->afbgchdm"(A, conj(A)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc1_leg3(FLo, ACu, M, ACd, FRo))
end

function contract_o1(FLo::leg3, ACu, A, ACd, FRo, O)
    D1,D2,D3,D4,_ = size(A)
    M = reshape(ein"(abcde,en),fghmn->afbgchdm"(A, O, conj(A)), D1^2,D2^2,D3^2,D4^2)
    return sum(oc1_leg3(FLo, ACu, M, ACd, FRo))
end

function contract_n1(FLo::leg4, ACu, A, ACd, FRo)
    return sum(oc1_leg4(FLo, ACu, A, conj(A), ACd, FRo))
end

function contract_o1(FLo::leg4, ACu, A, ACd, FRo, O)
    return sum(oc1_leg4(FLo, ACu, ein"abcde,ef->abcdf"(A, O), conj(A), ACd, FRo))
end

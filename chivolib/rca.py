
import numpy as np

def rca(inp, thresh, magic, on, off, centre, ):
    #allocate memory for out array
    out = np.empty_like(caa)

    shape = caa.shape

    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                """
                If the corresponding input pixel is equal to or greater than the
                magic value, copy it to the output.
                """
                if inp[i,j,k] >= magic:
                    out[i,j,k] = magic

                """
                If the corresponding input pixel is off, then the output must be also be
                off if "centre" is true.
                """
                if centre and inp[i,j,k] != on:
                    out[i,j,k] = off

                else:
                    sum = 0
                    tot = 0
                    for oi in range(i-1, i+2):
                        if oi < 0 or oi > shape[0] : continue
                        for oj in range(j-1, j+2):
                            if oj < 0 or oj > shape[1]: continue
                            for oz in range(z-1, z+2):
                                if oz <0 or oz > shape[2]: continue
                                tot += 1
                                if inp[oi,oj,oz] == on: sum += 1
                    if np.float(sum)/np.float(tot) > tresh :
                        out[i,j,k] = on
    return out

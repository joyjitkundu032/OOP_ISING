import numpy as np
import ising_oop as I

if __name__ == "__main__":

        import os
        try:
                os.mkdir("Data01")
        except OSError:
                pass

        for T in np.linspace(0.4,4.0,50):
		z=I.ising(20, T)
                z.initialize()
                z.cal_neighbour()
                z.cal_energy()
                f=open("Data01/Temp_%.4f.csv" %T,"a+")
                print(T)
                for i in np.arange(60000):
                        z.MC_steps()
                        M = z.mag
                        if(i % 200 == 0):
                                f.write("%d %d\n" % (i,M))
                f.close

import cpp_math
import cpp_rect
import matplotlib.pyplot as plt
import numpy as np

def fy420():
    return 420 #MPa

def fc30():
    return 30 #MPa

def sayHello():
    print("[Python] sayHello()!")

def add(a, b):
    print("[Python] defines function add(a,b)")
    return a + b

def calculate(n):
    print("[Python] call function sqr from C++")
    return cpp_math.multiplication(n)

def plots():
    print("plot from Python")
    plt.plot([1,2,3], [3,1,12])
    #return plt.show()

def point1():
    return 10.1

def point2():
    return 100.1

def get_area_rect(r):
    print("[Python] call class Rect from C++")
    return r.get_area()

#my_array = np.zeros(2)
def array__():
    my_array = np.array([[1,2,3], [11,22,33]])
    return my_array

def tuple_in_array():
    print("example list_of_tuples:")
    list_of_tuples = np.array([(101, 202, 303), (1001,2002,3003)])
    for (x,y,z) in list_of_tuples:
        print(x, y, z)
    print("end of example list_of_tuples")
    
    my_tup_arr = np.array([(100,200,300), (1100,2200,3300)])
    return my_tup_arr

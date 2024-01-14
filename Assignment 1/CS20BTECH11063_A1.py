def extendedEuclidean(a, b):
    A = a
    B = b
    x = 1
    y = 0
    u = 0
    v = 1
    
    while B != 0:
        q = A//B
        r = A%B
        tmp1 = x-u*q
        tmp2 = y-v*q
        
        A = B
        B = r
        x = u
        y = v
        u = tmp1
        v = tmp2
        
    return (x,y,A)

if __name__ == "__main__":
    # read input from file
    f = open("input-gcd.csv", "r")
    
    # read each line till end of file is reached
    for line in f:
        a, b = line.split(",")
        a = int(a)
        b = int(b)
        
        # compute gcd
        x, y, gcd = extendedEuclidean(a, b)
        
        # print the result
        # print("x = " + str(x) + ", y = " + str(y) + ", c = " + str(gcd))
        print("x = {}, y = {}, c = {}".format(x, y, gcd))
    
    f.close()
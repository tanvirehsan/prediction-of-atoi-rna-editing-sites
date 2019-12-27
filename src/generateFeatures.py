import itertools
import numpy as np

 
RNAelements = 'ACGU'
 

def sequenceType(seqType):
    
    if seqType == 'RNA':
        elements = RNAelements
    else:
        elements = None

    return elements


trackingFeatures = []

def genFeatures(pkTuple,psequenceType, X, Y):

    elements = sequenceType(psequenceType.upper())

    m2 = list(itertools.product(elements, repeat=2))
    m3 = list(itertools.product(elements, repeat=3))
    m4 = list(itertools.product(elements, repeat=4))
    m5 = list(itertools.product(elements, repeat=5))

    T = []  # All instance ...
    
    def countSeqFrequency(x):
        #4 Features        
        U = x.count('U')
        A = x.count('A');
        C = x.count('C');
        G = x.count('G');
        
        t.append(A)  
        t.append(C)  
        t.append(G)
        t.append(U)

    def kmers(seq, k):
        v = []
        for i in range(len(seq) - k + 1):
            v.append(seq[i:i + k])
        return v

    def pseudoKNC(x, k):
        ### k-mer ###
        ### A, AA, AAA

        for i in range(1, k + 1, 1):
            v = list(itertools.product(elements, repeat=i))
            # seqLength = len(x) - i + 1
            for i in v:
                # print(x.count(''.join(i)), end=',')
                t.append(x.count(''.join(i)))
        ### --- ###

    def zCurve(x):
        ### Z-Curve ### total = 3

            TU = x.count('U')

            A = x.count('A'); C = x.count('C'); G = x.count('G');

            x_ = (A + G) - (C + TU)
            y_ = (A + C) - (G + TU)
            z_ = (A + TU) - (C + G)
            # print(x_, end=','); print(y_, end=','); print(z_, end=',')
            t.append(x_); t.append(y_); t.append(z_)
            ### print('{},{},{}'.format(x_, y_, z_), end=',')
            ### --- ###
            # trackingFeatures.append('x_axis'); trackingFeatures.append('y_axis'); trackingFeatures.append('z_axis')

    def gcContent(x): 

            TU = x.count('U') 
            A = x.count('A');
            C = x.count('C');
            G = x.count('G');

            t.append( (G + C) / (A + C + G + TU)  * 100.0 )


    def cumulativeSkew(x):

            TU = x.count('U')
            A = x.count('A');
            C = x.count('C');
            G = x.count('G');

            GCSkew = (G-C)/(G+C)
            ATSkew = (A-TU)/(A+TU)

            t.append(GCSkew)
            t.append(ATSkew)


    def atgcRatio(x):

            TU = x.count('U') 
            A = x.count('A');
            C = x.count('C');
            G = x.count('G');

            t.append( (A+TU)/(G+C) )


    

    def generateFeatures(kTuple, x, y):
        
        countSeqFrequency(x) 
        zCurve(x)  
        gcContent(x) 
        cumulativeSkew(x)
        atgcRatio(x)
        pseudoKNC(x, kTuple)            #k=2|(16), k=3|(64), k=4|(256), k=5|(1024);
        # Features Generation Complete     
        ##############################################################
        #Target Variable
        t.append(y)
        #######################


    for x, y in zip(X, Y):
        t = []
        generateFeatures(pkTuple, x, y)
        T.append(t)
        ### --- ###

    ############################

    return np.array(T)



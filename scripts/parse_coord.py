import sys
import re

def main():
    infile = sys.argv[1]

    with open(infile) as f:
        lines = f.read().splitlines()
    # first we have to doctor the coords string
    # examples: 
    # nt:U50746.1   coords:4..867 has to become U50746.1:4..867
    # nt:U43883.1     join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,209..259)
    # has to become join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,U43883.1:209..259)

        for line in lines:
            try:
                aa = line.split('\t')
                acc = aa[0]
                coords = aa[1] 
                #print (acc, coords)


                identexpression = re.compile("[a-zA-Z0-9_.]+\:")
                startexpression = re.compile("[<>\d]+\.\.")
                stopexpression = re.compile("\.\.[<>\d]\d*")
        
                if coords.find("complement") != -1:
                    if coords.find("join") != -1: 
                        coords = coords.split(",")
                        for k in coords:
                            print (k)
                            ident = identexpression.search(k)
                            name = ident.group(0)
                            print (k, ident, name)
            except:
                pass
            #    print ("ERROR ! %s"%line )
            #    sys.exit()

    f.close()
   
if __name__ == '__main__':
    main()
    


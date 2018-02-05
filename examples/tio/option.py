"""
   Shows the action of the getOption method

Author: Will Hossack: The University of Edinburgh
"""
import tio as t

def zz(o,n):
    return n,o[n]

def main():

    opts = "Exit","Quit","Continue","Restart","Reset"

    ans = t.getOption("Option",opts)
    t.tprint("Option Number ",ans[0]," Option name ", ans[1])

    c = zz(opts,2)
    t.tprint(c)

main()

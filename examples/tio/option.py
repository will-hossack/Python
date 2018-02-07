"""
   Shows the action of the getOption method

Author: Will Hossack: The University of Edinburgh
"""
import tio as t

def main():

    while True: 
        opts = "exit","quit","continue","restart","reset"
        ans = t.getOption("Option",opts)
        t.tprint("Option Number ",ans[0]," Option name ", ans[1])
        if ans[0] == 0 or ans[0] == 1:
            break

main()

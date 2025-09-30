# -*- coding: utf-8 -*-
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser( prog="JANA-tools", description="commandline interface for JANA-tools")
    
    subparsers = parser.add_subparsers( title="subcommands", help="available tools" )
    
    ### RUNTIME analysis
    parsers_runtime = subparsers.add_parser("runtime",  help="runtime analysis")
    parsers_runtime.add_argument("-p", "--path", type=Path,  help="path to file or folder")
    parsers_runtime.add_argument("-l", "--logfile", default=False)   
    parsers_runtime.set_defaults(tool="runtime")
    
    ### BATCH
    parsers_batch = subparsers.add_parser("batch", help="batch modify and run dynamical refinements")
    #parsers_batch.add_argument("parameters", default=dict())
    parsers_batch.add_argument("-r", "--run_jana", action="store_true", help="Call and run JANA. Otherwise, instructions are given how to call Jana.")
    parsers_batch.add_argument("parameters", default=dict())
    parsers_batch.set_defaults(tool="batch")

        
    args = parser.parse_args()
    
    # TESTING / DEBUGGING
    # print(args)
    
    if not hasattr(args, "tool"):
        print("No tool selected.")
        # print more info ? help ?
    
    elif args.tool == 'runtime':
        from jana_tools.io_utils import RefFileHandler
        
        if args.path.is_file():
            # single ref file
            ref = RefFileHandler(args.path.with_suffix(".ref"))
            ref.report_runtime()
        elif args.path.is_dir():
            # search for REF files
            files = list( args.path.glob("*.ref") )
            
            times = []
            for file in files:
                ref = RefFileHandler(file)
                ref.report_runtime()
                times.append( ref.runtime_per_cycle )
            
            # Comprehensive listing at the end
            print("Average refinement runtime per cycle")
            print(f"{'File':<60}{'t (h:mm:ss)':>12}{'t/seconds':>12}")
            for file,t in zip(files,times):
                if t:
                    print(f"{file.stem:<60}{ref.duration_format(t):>12}{t:>12.1f}") 
                    
    #elif args.tool == 'batch':        
    else: 
        print("Not implemented yet.")


# testing, calling run.py directly
if __name__ == "__main__":
    main()


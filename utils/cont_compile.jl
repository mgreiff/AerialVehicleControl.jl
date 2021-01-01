function recompile()
    curpath = pwd();
    cd(CONT_LIB_SRC)
    run(`make clean`)
    run(`make run CFLAGS='-fPIC -shared -DTRACE_LEVEL=2'`)
    cd(curpath)
end

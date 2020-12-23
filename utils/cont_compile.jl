function recompile()
    curpath = pwd();
    cd(CONT_LIB_SRC)
    run(`make clean`)
    run(`make tests CFLAGS='-fPIC -shared'`)
    cd(curpath)
end

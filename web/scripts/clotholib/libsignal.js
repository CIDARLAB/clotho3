var libsignal = new Object();
libsignal._dispatch = new Object();

/* provide external implementation */
libtransport.emitSignal = function (signal, args) {
    if (signal in libsignal._dispatch) {
        libsignal._dispatch[signal](args);
    } else {
        alert("Bad signal: " + signal);
    }
}

libsignal._dispatch.socketIDChange = function (args) {
    /* TODO: is this needed? */
}

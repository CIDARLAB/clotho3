/* Needs to match up with
 *     `org.clothocad.core.aspects.Communicator.CLRegistry`
 */
var libsend = new Object();

libsend._init = function () {
    var x = null;
    var d = {
        "getTabList":x,
        "linkPage":x,
        "linkWidget":x,
        "serverEval":x,
        "submitCommand":x,
        "unlinkPage":x,
        "unlinkWidget":x,
        "updateQuery":x,
    };
    libsend._dispatch = d;
};

libsend._init();

libsend.call = function (channel, message) {
//alert("libsend calling "+message+" in channel "+channel);
    if (channel in libsend._dispatch) {
        var auth_key = "";
        try {
            auth_key = libcookie.get("auth_key");
        } catch (Error) {}
        var flags = {"auth_key": auth_key};
        libtransport.send(channel, message, JSON.stringify(flags));
    } else {
        alert("Cannot send message: channel " + channel + " is not defined");
    }
};

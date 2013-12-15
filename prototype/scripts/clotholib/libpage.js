/* TODO: register BROWSE and EDITOR modes */
var libpage = new Object();

/* need external implementation */
libpage.getCurrentMode = function () {
    alert("Missing external implementation");
}

/* Associates page modes with URLs
 * Reflects org.clothocad.core.aspects.Communicator.Mind.PageMode on server
 */
libpage.modes = {
    "BIO"    : "/servlet/page?BIO",
    "BROWSE"    : "/servlet/page?BROWSE",
    "EDITOR"    : "/servlet/page?EDITOR",
    "HOMEPAGE"  : "/servlet/page?HOMEPAGE",
    "MYSTUFF"   : "/servlet/page?MYSTUFF",
    "TRAILS"    : "/servlet/page?TRAILS",
    "WORKSPACE" : "/servlet/page?WORKSPACE"
};

libpage.add = function (mode, ephemeral_link_page_id) {
    var window_open_args = ["toolbar=1",
                            "location=1",
                            "directories=1",
                            "toolbar=1",
                            "status=1",
                            "menubar=1",
                            "scrollbars=1",
                            "resizeable=1",
                           ].join(",");
    window.open(libpage.modes[mode], ephemeral_link_page_id, window_open_args);
};

libpage.gofullscreen = function () {};

libpage.getdivbyuuid = function (uuid) { return $("[uuid="+uuid+"]");};

/* Widget API */

var API = new Object();
API.widget = new Object();

API.show = function (view_uuid, args) {
    libsend.call("serverEval",
                 "clotho.show('" + view_uuid + "', '" + args + "');");
}

API.widget.getParent = function (domobj) {
    do {
        domobj = domobj.parentElement;
    } while (domobj !== null && !domobj.hasAttribute("widget_id"));
    return (domobj === null) ? null : domobj;
}

API.widget.getID = function (domobj) {
    return (domobj === null || !domobj.hasAttribute("widget_id")) ?
           null : domobj.getAttribute("widget_id");
}

API.widget.getByID = function (id) {
    var nodelist = document.getElementsByTagName("div");
    for (var i = 0; i < nodelist.length; i++) {
        if (nodelist[i].getAttribute("widget_id") == id) return nodelist[i];
    }
    return null;
}

API.widget.call = function (domobj, methname, args) {
    domobj[methname](domobj, args);
}

API.widget.parentCall = function (domobj, methname, args) {
    API.widget.call(API.widget.getParent(domobj), methname, args);
}

API.widget.getContainer = function (domobj) {
    for (var child = domobj.firstChild;
         child !== null && child !== undefined;
         child = child.nextSibling) {
        if (child.getAttribute !== undefined &&
            child.getAttribute("class") === "widget-container") {
            return child;
        }
    }
}

/* Derived from: http://www.w3schools.com/js/js_cookies.asp
 *               http://techpatterns.com/downloads/javascript_cookies.php
 */
var libcookie = new Object();

libcookie.get = function (c_name) {
    var cookies = document.cookie.split(";");
    for (i in cookies) {
        var c = cookies[i].split("=");
        if (c.length != 2) { throw Error(); }
        if (c[0].replace(/^\s+|\s+$/g, "") == c_name) {
            return unescape(c[1]);
        }
    }
    throw Error();
}

libcookie.set = function (c_name, value) {
    var c_value = escape(value) + "; path='/'";
    document.cookie = c_name + "=" + c_value;
}

libcookie.del = function (c_name) {
    document.cookie = c_name + "; expires=01-Jan-1970 00:00:00";
}

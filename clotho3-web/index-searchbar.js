'use strict';

/*!
 Async script loader
 $script.js v1.3
 usage: https://github.com/ded/script.js
 */
!function(a,b,c){function t(a,c){var e=b.createElement("script"),f=j;e.onload=e.onerror=e[o]=function(){e[m]&&!/^c|loade/.test(e[m])||f||(e.onload=e[o]=null,f=1,c())},e.async=1,e.src=a,d.insertBefore(e,d.firstChild)}function q(a,b){p(a,function(a){return!b(a)})}var d=b.getElementsByTagName("head")[0],e={},f={},g={},h={},i="string",j=!1,k="push",l="DOMContentLoaded",m="readyState",n="addEventListener",o="onreadystatechange",p=function(a,b){for(var c=0,d=a.length;c<d;++c)if(!b(a[c]))return j;return 1};!b[m]&&b[n]&&(b[n](l,function r(){b.removeEventListener(l,r,j),b[m]="complete"},j),b[m]="loading");var s=function(a,b,d){function o(){if(!--m){e[l]=1,j&&j();for(var a in g)p(a.split("|"),n)&&!q(g[a],n)&&(g[a]=[])}}function n(a){return a.call?a():e[a]}a=a[k]?a:[a];var i=b&&b.call,j=i?b:d,l=i?a.join(""):b,m=a.length;c(function(){q(a,function(a){h[a]?(l&&(f[l]=1),o()):(h[a]=1,l&&(f[l]=1),t(s.path?s.path+a+".js":a,o))})},0);return s};s.get=t,s.ready=function(a,b,c){a=a[k]?a:[a];var d=[];!q(a,function(a){e[a]||d[k](a)})&&p(a,function(a){return e[a]})?b():!function(a){g[a]=g[a]||[],g[a][k](b),c&&c(d)}(a.join("|"));return s};var u=a.$script;s.noConflict=function(){a.$script=u;return this},typeof module!="undefined"&&module.exports?module.exports=s:a.$script=s}(this,document,setTimeout);

//future - use CDN for these
$script("lib/jquery-1.9.1.js", "jquery", function() {
    $script("lib/jquery-ui-1.10.2.js", "jquery-ui");
});

//note - dependencies prevent full asynchrony
//todo - better jquery handling... just need to be ready before DOM content ready
$script.ready("jquery", function() {
    $script("lib/angular-1.1.5/angular.js", "angular", function() {

        var Application = window.Application = {};
        Application.Search = angular.module('clotho.search', ['clotho.core', 'clotho.interface']);

        Application.Foundation = angular.module('clotho.core', [])
            .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {
                //on first run, add API to $clotho object
                window.$clotho.api = Clotho;

                //extend scope with Clotho API
                $rootScope.Clotho = Clotho;

                $rootScope.$safeApply = function(fn) {
                    fn = fn || function() {};
                    if($rootScope.$$phase) { fn(); }
                    else { $rootScope.$apply(fn); }
                };
            }]);

        //dummy function, won't be used in this context
        Application.Interface = angular.module('clotho.interface', ['clotho.core']);
        Application.Interface.provider("$dialog", function(){
            this.$get = ['$rootScope', function() {
                return {}
            }];
        });

        //load all angular dependencies asynchronously
        $script([
            //foundation
            "foundation/clientAPI.js",
            "foundation/clothoAPI.js",
            "foundation/pubsub-service.js",
            "foundation/socket-service.js",
            "foundation/collector-service.js",

            //UI
            "interface/ui-services.js",
            "interface/ui-directives.js",


            //search
            "search/search-service.js",
            "search/search-directives.js"
        ], "angular-stack");
    });
});

$script.ready(["jquery", "angular-stack"], function() {
    var $clotho = window.$clotho = {};

    $clotho.socket = new WebSocket("ws://localhost:8080/websocket");
    $clotho.socket.onopen = function () {
        angular.bootstrap(
            angular.element(document),
            ['clotho.search']);
    }
});
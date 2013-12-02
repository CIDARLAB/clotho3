angular.module("clotho.foundation",["clotho.core","clotho.setup","clotho.clothoDirectives","clotho.extensions"]),angular.module("clotho.setup",[]).run(["$rootScope","Clotho",function(a,b){a.Clotho=b}]),angular.module("clotho.clothoDirectives",[]).directive("clothoRun",["Clotho",function(a){var b={input:!0,textarea:!0,select:!0};return{restrict:"A",require:"ngModel",scope:!0,link:function(c,d,e,f){var g=!!b[angular.lowercase(d[0].nodeName)];g&&(f.$render=angular.noop);var h=!1;c.$watch(function(){return e.clothoRun},function(a){a&&l(f.$modelValue)}),c.$watch(function(){return f.$modelValue},function(a){l(a)}),c.$watch(function(){return e.clothoRunUpdateModel},function(a){h=!!a});var i=function(a){return angular.isArray(a)?a:[a]},j=function(a){h&&f.$setViewValue(a)},k=function(a){var b=g?"val":"text";d[b](a)},l=function(b){return b=i(b),a.run(e.clothoRun,b).then(function(a){console.log(a),j(a),k(a)})}}}}]),!function(a,b,c){function d(a,c){var d=b.createElement("script"),e=k;d.onload=d.onerror=d[p]=function(){d[n]&&!/^c|loade/.test(d[n])||e||(d.onload=d[p]=null,e=1,c())},d.async=1,d.src=a,f.insertBefore(d,f.firstChild)}function e(a,b){q(a,function(a){return!b(a)})}var f=b.getElementsByTagName("head")[0],g={},h={},i={},j={},k=!1,l="push",m="DOMContentLoaded",n="readyState",o="addEventListener",p="onreadystatechange",q=function(a,b){for(var c=0,d=a.length;d>c;++c)if(!b(a[c]))return k;return 1};!b[n]&&b[o]&&(b[o](m,function t(){b.removeEventListener(m,t,k),b[n]="complete"},k),b[n]="loading");var r=function(a,b,f){function k(){if(!--s){g[p]=1,o&&o();for(var a in i)q(a.split("|"),m)&&!e(i[a],m)&&(i[a]=[])}}function m(a){return a.call?a():g[a]}a=a[l]?a:[a];var n=b&&b.call,o=n?b:f,p=n?a.join(""):b,s=a.length;return c(function(){e(a,function(a){j[a]?(p&&(h[p]=1),k()):(j[a]=1,p&&(h[p]=1),d(r.path?r.path+a+".js":a,k))})},0),r};r.get=d,r.ready=function(a,b,c){a=a[l]?a:[a];var d=[];return!e(a,function(a){g[a]||d[l](a)})&&q(a,function(a){return g[a]})?b():!function(a){i[a]=i[a]||[],i[a][l](b),c&&c(d)}(a.join("|")),r};var s=a.$script;r.noConflict=function(){return a.$script=s,this},"undefined"!=typeof module&&module.exports?module.exports=r:a.$script=r}(this,document,setTimeout),angular.module("clotho.extensions",[]).config(["$routeProvider","$controllerProvider","$compileProvider","$filterProvider","$provide",function(a,b,c,d,e){window.$clotho.extensions=$clotho.extensions={},$clotho.extensions.providers={$routeProvider:a,$controllerProvider:b,$compileProvider:c,$filterProvider:d,$provide:e},$clotho.extensions.getQueue=function(){return angular.module("clotho.extensions")._invokeQueue},$clotho.extensions.registeredQueue=$clotho.extensions.getQueue().length}]).run(["$rootScope","$q","$timeout","$templateCache","$http","$rootElement",function(a,b,c,d,e,f){$clotho.extensions.processQueue=function(){for(var a=$clotho.extensions.getQueue(),b=$clotho.extensions.registeredQueue;b<a.length;b++){var c=a[b],d=$clotho.extensions.providers[c[0]];d&&d[c[1]].apply(d,c[2])}$clotho.extensions.registeredQueue=$clotho.extensions.getQueue().length},$clotho.extensions.recompile=function(a,b){"undefined"!=typeof a&&(b=b||{},f.injector().invoke(function(c,d){var e=d.$new();angular.extend(e,b),c($(a))(e),d.$apply()}))},$clotho.extensions.mixin=function(c,d,e){if(angular.isUndefined(c)||""==c)return b.when("no mixin url");var f=b.defer();return $script(c,function(){$clotho.extensions.processQueue(),$clotho.extensions.recompile(d,e),a.$safeApply(f.resolve(d))}),f.promise},$clotho.extensions.script=function(c){if(angular.isUndefined(c)||""==c)return b.when("no script url");var d,e=b.defer();return angular.isString(c)?d=c+"?_="+Date.now():(d=[],angular.forEach(c,function(a){d.push(a+"?_="+Date.now())})),$script(c,function(){a.$safeApply(e.resolve())}),e.promise},$clotho.extensions.css=function(c){if(angular.isUndefined(c)||""==c)return b.when("no css url");var d=b.defer();if(document.createStyleSheet)document.createStyleSheet(c),a.$safeApply(d.resolve());else{var e=document.createElement("link");e.type="text/css",e.rel="stylesheet",e.href=c,document.getElementsByTagName("head")[0].appendChild(e),a.$safeApply(d.resolve())}return d.promise},$clotho.extensions.cache=function(a){if(angular.isUndefined(a)||""==a)return b.when();var c=b.defer();return e.get(a,{cache:d}).success(function(){c.resolve()}).error(function(){c.reject()}),c.promise};var g=0;$clotho.extensions.bootstrap=function(a){g++;var c=b.defer(),d=$($('<div clotho-widget clotho-widget-uuid="'+g+'" clotho-widget-name="'+a.moduleName+'"></div>').append("<div ng-view></div>")).appendTo($clotho.appWidgets);return $clotho.extensions.script(a.moduleUrl).then(function(){angular.bootstrap(d,[a.moduleName]),c.resolve([g,"[clotho-widget-uuid="+g+"]"])}),c.promise}}]);
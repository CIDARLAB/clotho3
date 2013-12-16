var $clotho=window.$clotho=window.$clotho||{};angular.module("clotho.core",["clotho.angularAdditions"]),angular.module("clotho.angularAdditions",["ngSanitize","ngRoute"]).config(function(){var a={};a.isEmpty=function(a){return angular.isUndefined(a)||""===a||null===a||a!==a},a.isScope=function(a){return a&&a.$evalAsync&&a.$watch},angular.extend(angular,a)}).run(["$rootScope",function(a){a.$safeApply=function(){var a,b,c=!1;if(1==arguments.length){var d=arguments[0];"function"==typeof d?b=d:a=d}else a=arguments[0],b=arguments[1],3==arguments.length&&(c=!!arguments[2]);a=a||this,b=b||function(){},c||!a.$$phase?a.$apply?a.$apply(b):a.apply(b):b()}}]),angular.module("clotho.core").service("Clotho",["Socket","Collector","PubSub","$window","$q","$rootScope","$location","$timeout",function(a,b,c,d,e,f,g,h){function i(){var d={};d.send=function(b){a.send(angular.toJson(b))},d.pack=function(a,b,c,d){return{channel:a,data:b,requestId:c,options:d}},d.emit=function(a,b,c,e){d.send(d.pack(a,b,c,e))},d.api={},d.searchbar={},d.api.emit=function(a,b,c){d.emit(a,b,c)},d.searchbar.emit=function(a,b,c){d.emit(a,b,c)},d.emitSubCallback=function(a,b,g,i){var j=e.defer(),k=Date.now().toString();angular.isFunction(g)||(g=angular.noop());var l=h(function(){j.reject(null)},5e3);return c.once(a+":"+k,function(a){h.cancel(l),f.$safeApply(j.resolve(a)),g(a)},"$clotho"),d.emit(a,b,k,i),j.promise},d.emitSubOnce=function(a,b,c){return d.emitSubCallback(a,b,angular.noop,c)};var i=function(a,b){var c={username:a,password:b};return d.emitSubOnce("login",c)},j=function(a,b,c){return c="undefined"!=typeof c?!!c:!1,c?l(a,b):k(a,b)},k=function(a,c){var f=b.retrieveModel(a);if(f){var g=e.defer();return g.resolve(f),g.promise}var h=function(c){b.storeModel(a,c)};return d.emitSubCallback("get",a,h,c)},l=function(a){var c=b.retrieveModel(a);return c?c:void 0},m=function(a){if(angular.isEmpty(a))return!1;for(;a.$$v;)a=a.$$v;var c=function(){b.storeModel(a.id,a)};return d.emitSubCallback("set",a,c)},n=function(a,b,d){d="undefined"!=typeof d?d:null,c.on("update:"+a,function(a){f.$safeApply(b(a))},d)},o=function(a,b,d,e){e="undefined"!=typeof e?e:null,c.on("update:"+a,function(a){f.$safeApply(b[d]=a)},e)},p=function(a,b,d){d="undefined"!=typeof d?d:null,c.on(a,function(a){f.$safeApply(b(a))},d)},q=function(a,b){c.trigger(a,b)},r=function(a){c.destroy(a)},s=function(a,b){d.api.emit(a,b||{})},t=function(a,b){d.api.emit("broadcast",d.pack(a,b))},u=function(b,c){return a.on(b,c)},v=function(b,c){return a.once(b,c)},w=function(b){a.off(b)},x=function(a,c){var e=function(a){angular.forEach(a,function(a){b.storeModel(a.id,a)})};return d.emitSubCallback("query",a,e,c)},y=function(a){return d.emitSubOnce("create",a)},z=function(a){var c=function(){b.removeModel(a)};return d.emitSubCallback("destroy",a,c)},A=function(a){g.path("/editor/"+a)},B=function(a,b){var c={uuid:a,timestamp:b};return d.emitSubOnce("revert",c)},C=function(a){return d.emitSubOnce("validate",a)},D=function(a,b){var c={uuid:a,args:b};d.api.emit("show_old",c)},E=function(a){d.api.emit("show",a)},F=function(){console.log("need to set up share")},G=0,H=function(a){G++;var b=G,c=e.defer(),d=$($('<div clotho-widget clotho-widget-uuid="'+b+'" clotho-widget-name="'+a.moduleName+'"></div>').append("<div ng-view></div>")).appendTo($clotho.appWidgets),g=a.moduleUrl?a.moduleUrl+"?_="+Date.now():"";return $script(g,function(){angular.bootstrap(d,[a.moduleName]),c.resolve([b,"[clotho-widget-uuid="+b+"]"]),f.$safeApply()}),c.promise},I=function(a){d.api.emit("log",a)},J=function(a,b){var c={userID:b||"",msg:a,timestamp:Date.now()};d.api.emit("say",c)},K=function(a){d.api.emit("notify",a)},L=function(a,b){var c={userID:b,msg:a};d.api.emit("alert",c)},M=function(a){var b={query:a};return d.emitSubOnce("autocomplete",b)},N=function(a){var g=e.defer();if(b.hasItem("detail_"+a))g.resolve(b.retrieveModel("detail_"+a));else{var h={uuid:a};d.searchbar.emit("autocompleteDetail",h),c.once("update:detail_"+a,function(a){f.$safeApply(g.resolve(a))},"$clotho")}return g.promise},O=function(a){return d.emitSubOnce("submit",a)},P=function(a,b){var c={id:a,args:b};return d.emitSubOnce("run",c)},Q=function(a){return d.emitSubOnce("recent",a)},R=function(a,b,c){return P("gradeQuiz",[a,b,c])};return{login:i,get:j,set:m,query:x,create:y,edit:A,validate:C,revert:B,destroy:z,show:E,show_old:D,say:J,log:I,alert:L,run:P,recent:Q,notify:K,gradeQuiz:R,bootstrap:H,watch:n,watch2:o,listen:p,silence:r,trigger:q,emit:s,broadcast:t,on:u,once:v,off:w,share:F,submit:O,autocomplete:M,autocompleteDetail:N}}return d.$clotho.api?d.$clotho.api:d.$clotho.api=i()}]),angular.module("clotho.core").service("ClientAPI",["PubSub","Collector","$q","$templateCache","$http","$rootScope","$location","$compile",function(a,b,c,d,e,f,g){var h,i=!1;try{angular.module("clotho.interface"),i=!0,h=angular.injector("clotho.interface").get("$dialog")}catch(j){}var k=function(a){var c=a,d=c.id||c.uuid||!1;b.storeModel(d,c)},l=function(a){if(i){var b={backdrop:!0,keyboard:!0,backdropClick:!0,template:'<form sharable-editor name="sharableEditor" uuid="'+a+'" class="span6 form-horizontal well" novalidate></form>'},c=h.dialog(b);c.open()}else g.path("/editor/"+a)},m=function(a){g.path(a)},n=function(){},o=function(b){a.trigger(b.channel,b.data)},p=function(a,b){b.model,b.template},q=function(a){console.log(a);var b=a.template,c=a.controller||"",g=a.args||{},h=a.dependencies||[],i=a.styles||{},j=a.target&&$($clotho.appRoot).has(a.target)?a.target:$($clotho.appRoot).find("[ng-view]");f.$safeApply(e.get(b,{cache:d}).success(function(a){$clotho.extensions.mixin([h,c],$(a).appendTo(j),g).then(function(a){a.css(i)})}).error(function(){console.log("error getting template")}))},r=function(){},s=function(a){console.log(a)},t=function(b){b.timestamp=b.timestamp?b.timestamp:Date.now(),b.from=b.from?b.from:"server",a.trigger("activityLog",b)},u=function(b){a.trigger("serverAlert",{}),i?f.$safeApply(h.serverAlert(b).open().then(function(a){console.log("dialog closed with result: "+a)})):window.alert(b)},v=function(){},w=function(b,c){a.trigger("revisions:"+b,c)},x=function(a){g.path("/trails/"+a)},y=function(c){console.log("Hit autocompletedetail"),console.log(c);var d=c.command_object.function_id;b.storeModel("detail_"+d,c),a.trigger("autocompleteDetail_"+d,c)};return{collect:k,edit:l,changeUrl:m,notify:n,broadcast:o,log:s,say:t,alert:u,display:q,display_simple:q,display_old:p,hide:r,help:v,revisions:w,startTrail:x,autocompleteDetail:y}}]),angular.module("clotho.core").service("Collector",["$window","$document","PubSub",function(a,b,c){function d(){function b(a,b){c.trigger("update",a),c.trigger("update:"+a,angular.copy(b))}function d(){var c=a.localStorage,d=JSON,e="clotho_",f=function(){try{return"localStorage"in a&&null!==a.localStorage}catch(b){return!1}},g=function(){return c.clear(),this},h=function(a,b){var f=c.getItem(e+a);return"undefined"==typeof f||f===!1?"undefined"!=typeof b?b:null:d.parse(f)},i=function(a){var b=c.getItem(e+a);return null!=b&&"undefined"!=typeof b},j=function(a){return c.removeItem(e+a),this},k=function(a,b){return b||0===b||""===b?(c.setItem(e+a,d.stringify(b)),this):!1},l=function(c){c||(c=a.event);var d=c.key.replace(e,"")||"",f=h(d);console.log("handle_storage_event for "+d),b(d,f)};return a.addEventListener?a.addEventListener("storage",l,!1):a.attachEvent("onstorage",l),{getPrefix:function(){return e},isSupported:f,clear:g,getItem:h,hasItem:i,removeItem:j,setItem:k}}var e=d(),f={},g=function(a,b){f[a]=b,e.setItem(a,b)},h=function(a,c,d){d||!angular.equals(f[a],c)?(g(a,c),b(a,c)):console.log("COLLECTOR	"+a+"model is same as collector")},i=function(a){return angular.copy(j(a))},j=function(a){return"undefined"==typeof a?f:f[a]&&"undefined"!=typeof f[a]?f[a]:(f[a]=e.getItem(a))?f[a]:!1},k=function(a){return f[a]?(f[a]=null,e.removeItem(a),!0):!1},l=function(){e.clear(),f={},c.trigger("collector_reset",null)};return{collector:f,hasItem:e.hasItem,silentAddModel:g,storeModel:h,retrieveModel:i,retrieveRef:j,removeModel:k,clearStorage:l}}return a.$clotho.$collector?a.$clotho.$collector:a.$clotho.$collector=d()}]),angular.module("clotho.core").service("PubSub",function(){function a(){var a={},b={},c=/\s+/,d=function(a,c){b[a]||(b[a]=[]),b[a].push(c)},e=function(b,c,e,f){a[b]||(a[b]=[]),a[b].push([c,f]);var g=[b,c];return e&&"undefined"!=typeof e&&(angular.isScope(e)?(e.$on("$destroy",function(){l(e.$id)}),d(e.$id,g)):d(e,g)),g},f=function(b){var c=b[0];a[c]&&angular.forEach(a[c],function(d,e){d[0]==b[1]&&a[c].splice(e,1)})},g=function(){console.log(a)},h=function(b,d){var e=b.split(c);angular.forEach(e,function(b){a[b]&&angular.forEach(a[b],function(c,e){c[0](d),c[1]&&a[b].splice(e,1)})})},i=function(a,b,d,f){if(d=d||null,f=f||!1,c.test(a)){var g=[],h=a.split(c);return angular.forEach(h,function(a){var c=e(a,b,d,f);g.push(c)}),g}return e(a,b,d,f)},j=function(a,b,c){i(a,b,c,!0)},k=function(a){angular.isArray(a[0])?angular.forEach(a,function(a){f(a)}):f(a)},l=function(a){var c=b[a];angular.forEach(c,function(a){f(a)}),b[a]=[]},m=function(d){"all"==d&&(a={},b={});var e=d.split(c);angular.forEach(e,function(b){a[b]=[]})};return{logListeners:g,trigger:h,on:i,once:j,off:k,destroy:l,clear:m}}return window.$clotho.$pubsub?window.$clotho.$pubsub:window.$clotho.$pubsub=a()}),angular.module("clotho.core").service("Socket",["$window","$q","PubSub","ClientAPI",function(a,b,c,d){function e(){var b,e,f=[];return b=a.$clotho.socket?a.$clotho.socket:a.$clotho.socket=new WebSocket("wss://"+window.location.host+window.location.pathname+"websocket"),1==b.readyState?(console.log("socket already present, sending items in socket Queue"),e=!0,angular.forEach(f,function(a,c){b.send(a),f.splice(c,1)})):b.onopen=function(){console.log("socket opened, sending queued items"),e=!0,angular.forEach(f,function(a,c){b.send(a),f.splice(c,1)})},b.onopen=function(){console.log("socket opened, sending queued items"),e=!0,angular.forEach(f,function(a){b.send(a)})},b.onerror=function(){e=!1,console.log("socket error")},b.onclose=function(){d.say({"class":"error",text:"Socket Connection Closed",from:"client"})},b.onmessage=function(a){a=JSON.parse(a.data),console.log("SOCKET	received",a);var b=a.channel,e=a.requestId,f=a.data;"function"==typeof d[b]?d[b](f):c.trigger(b+":"+e,f)},{emit:function(a,c,d){console.log("SOCKET	data emitted on channel: "+a);var e={channel:a,data:c};b.send(e),"function"==typeof d&&d(e)},send:function(a){return e?(console.log("SOCKET	sending data: "+a),b.send(a),void 0):(console.log("socket not ready, queueing request"),f.push(a),void 0)}}}return a.$clotho.$socket?a.$clotho.$socket:a.$clotho.$socket=e()}]);
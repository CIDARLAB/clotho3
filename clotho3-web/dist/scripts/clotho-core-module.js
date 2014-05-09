var $clotho=window.$clotho=window.$clotho||{};angular.module("clotho.core",["clotho.angularAdditions"]),angular.module("clotho.angularAdditions",[]).config(function(){var a={};a.isEmpty=function(a){if(angular.isNumber(a))return a===a;if(a===!0||a===!1)return!1;if(angular.isObject(a)){if(0===a.length)return!0;for(var b in a)if(a.hasOwnProperty(b))return!1;return!0}return angular.isUndefined(a)||null===a||0==a.length?!0:!1},a.isScope=function(a){return!!a&&angular.isFunction(a.$evalAsync)&&angular.isFunction(a.$watch)},a.once=function(a){var b,c;if(!angular.isFunction(a))throw new TypeError;return function(){return b?c:(b=!0,c=a.apply(this,arguments),a=null,c)}},a.remove=function(a,b,c){var d=-1,e=a?a.length:0,f=[];for(c=c||null;++d<e;){var g=a[d];(angular.isUndefined(g)||b.call(c,g,d,a))&&(angular.isDefined(g)&&f.push(g),a.splice(d--,1),e--)}return f},a.map=function(a,b,c){var d=angular.isArray(a)?[]:{};return angular.forEach(a,function(a,e,f){angular.isFunction(b)&&(angular.isArray(d)?d.push(b.call(c,a,e,f)):d[e]=b.call(c,a,e,f))}),d},angular.extend(angular,a)}).run(["$rootScope",function(a){a.$safeApply=function(){var a,b,c=!1;if(1==arguments.length){var d=arguments[0];"function"==typeof d?b=d:a=d}else a=arguments[0],b=arguments[1],3==arguments.length&&(c=!!arguments[2]);a=a||this,b=b||function(){},c||!a.$$phase?a.$apply?a.$apply(b):a.apply(b):b()}}]),angular.module("clotho.core").config(["$provide",function(a){a.decorator("$log",["$delegate",function(a){return a.table=function(){var b=[].slice.call(arguments);window.console&&window.console.table?console.table(b[0],b[1]):a.log(null,b)},a}])}]).factory("Debug",["$log","$window",function(a,b){function c(a){return d===!0?!1:e.indexOf(a)<0}var d=!1,e=[],f={};return function(d,e){d=angular.isString(d)?d:"ANONYMOUS",e=/(^#[0-9A-F]{6}$)|(^#[0-9A-F]{3}$)/i.test(e)?e:"#"+Math.floor(16777215*Math.random()).toString(16),f[d]=[];var g={};return angular.forEach(["log","warn","error","debug","info"],function(b){g[b]=function(){var g=new Date,h=g.getHours()+":"+g.getMinutes()+":"+g.getSeconds()+"."+g.getMilliseconds();f[d].push({message:arguments,time:g.valueOf()}),c(d)&&a[b].apply(null,["%c"+h+" - "+d+"	","color: "+e+";"].concat(Array.prototype.slice.call(arguments,0)))}}),angular.forEach(["table"],function(a){g[a]=b.console[a]}),g.object=function(b){f[d].push({message:b,time:Date.now().valueOf()}),angular.isObject(b)?a.log("%c"+JSON.stringify(b,null,2),"color: "+e+";"):a.log("%c"+b,"color: "+e+";")},g.$log=a,g.console=b.console||{},g}}]),angular.module("clotho.core").service("Clotho",["Socket","Collector","PubSub","Debug","$window","$q","$rootScope","$location","$timeout",function(a,b,c,d,e,f,g,h,i){function j(){var e=new d("ClothoAPI","#cc5555"),g={};g.send=function(b){a.send(angular.toJson(b))},g.pack=function(a,b,c,d){return{channel:a,data:b,requestId:c,options:d}},g.emit=function(a,b,c,d){g.send(g.pack(a,b,c,d))};var j=new function(){var a=0;this.next=function(){return a+=1,a.toString()}};g.emitSubCallback=function(a,b,d,e){var h=f.defer(),k=Date.now().toString()+j.next();angular.isFunction(d)||(d=angular.noop);var l=i(function(){h.reject(null)},5e3);return c.once(a+":"+k,function(a){i.cancel(l),angular.isUndefined(a)?h.reject(null):h.resolve(a),d(a)},"$clotho"),g.emit(a,b,k,e),h.promise},g.emitSubOnce=function(a,b,c){return g.emitSubCallback(a,b,angular.noop,c)};var k=function(a,b){var c={username:a,password:b};return g.emitSubOnce("login",c)},l=function(){return g.emitSubOnce("logout","")},m=function(a,b,c){return c=angular.isDefined(c)?!!c:!1,c?o(a,b):n(a,b)},n=function(a,b){return angular.isUndefined(a)?f.when():g.emitSubOnce("get",a,b)},o=function(a){if(!angular.isUndefined(a)){var c=b.retrieveModel(a);return c?c:void 0}},p=function(a){if(angular.isEmpty(a))return!1;for(;a.$$v;)a=a.$$v;var c=function(){b.storeModel(a.id,a)};return g.emitSubCallback("set",a,c)},q=function(a,b,d){return d="undefined"!=typeof d?d:null,c.on("update:"+a,function(a){angular.isFunction(b)?b.apply(d,[a]):angular.extend(b,a)},d)},r=function(a,b,d){d="undefined"!=typeof d?d:null,c.on(a,function(a){b(a)},d)},s=function(a,b){c.trigger(a,b)},t=function(a){c.destroy(a)},u=function(a,b){g.emit(a,b||{})},v=function(a,b){g.emit("broadcast",g.pack(a,b))},w=function(b,c){return a.on(b,c)},x=function(b,c){return a.once(b,c)},y=function(b){a.off(b)},z=function(a,c){var d=function(c){console.groupCollapsed("Query Results for: "+JSON.stringify(a));try{angular.forEach(c,function(a){b.storeModel(a.id,a)})}catch(d){e.warn("error saving all models",d)}console.groupEnd()};return g.emitSubCallback("query",a,d,c)},A=function(a){return g.emitSubOnce("create",a)},B=function(a){if(a){var c=function(){b.removeModel(a)};return g.emitSubCallback("destroy",a,c)}},C=function(a){h.path("/editor?id="+a)},D=function(a,b){var c={uuid:a,timestamp:b};return g.emitSubOnce("revert",c)},E=function(a){return g.emitSubOnce("validate",a)},F=function(a,b){var c={viewId:a,options:b};g.emit("show",c)},G=function(){console.log("need to set up share")},H=function(a){g.emit("log",a)},I=function(a,b){var c={userID:b||"",msg:a,timestamp:Date.now()};g.emit("say",c)},J=function(a){g.emit("notify",a)},K=function(a,b){var c={userID:b,msg:a};g.emit("alert",c)},L=function(a,b){var c={query:a};return g.emitSubOnce("autocomplete",c,b)},M=function(a){return angular.isString(a)&&(a={query:a,tokens:[]}),g.emitSubOnce("submit",a)},N=function(a,b,c){var d={id:a,args:b};return c=angular.isObject(c)?c:{},g.emitSubOnce("run",d,c)},O=function(a){return g.emitSubOnce("recent",a)};return{login:k,logout:l,get:m,set:p,query:z,create:A,edit:C,validate:E,revert:D,destroy:B,show:F,say:I,log:H,alert:K,run:N,recent:O,notify:J,submit:M,autocomplete:L,watch:q,listen:r,silence:t,trigger:s,emit:u,broadcast:v,on:w,once:x,off:y,share:G}}return e.$clotho.api?e.$clotho.api:e.$clotho.api=j()}]),angular.module("clotho.core").service("ClientAPI",["PubSub","Collector","Debug","$q","$injector","$window","$templateCache","$http","$rootScope","$location","$clothoModal","$document",function(a,b,c,d,e,f,g,h,i,j,k){var l=new c("ClientAPI","#dd99dd"),m=function(a){var c=a,d=c.id||c.uuid||!1;b.storeModel(d,c)},n=function(b){a.trigger(b.channel,b.data)},o=function(a,b){var c=angular.element('<div clotho-show="'+a+'"></div>'),d=document.querySelector(b);d||(d=document.getElementById("clothoAppWidgets")),d=angular.element(d),d.append($compile(c)(i))},p=function(a,b){var c=angular.element(document.querySelector('[clotho-show="'+a+'"]'));b.apply(null,c.remove())},q=function(a){l.info(a)},r=function(b){angular.isString(b.text)||(angular.isNumber(b.text)?b.text=parseInt(b.text):null===b.text?b.text="null":angular.isUndefined(b.text)?b.text="undefined":b.text===!0?b.text="true":b.text===!1&&(b.text="false"));var c={"class":"muted",from:"server",timestamp:Date.now().valueOf()};b=angular.extend({},c,b);var d={success:"success",warning:"warning",error:"danger",failure:"danger",normal:"success",muted:"muted",info:"info"};b.class=d[angular.lowercase(b.class)],a.trigger("activityLog",b)},s=function(){k.create({title:"Clotho Alert",content:"msg"})},t=function(){},u=function(b,c){a.trigger("revisions:"+b,c)},v=function(a){j.path(a)},w=function(a){e.has("$route")?j.path("/editor?id="+a):r({text:"Editor not available",from:"client","class":"error"})},x=function(a){j.path("/trails/"+a)};return{collect:m,broadcast:n,log:q,say:r,alert:s,display:o,hide:p,help:t,revisions:u,changeUrl:v,edit:w,startTrail:x}}]),angular.module("clotho.core").service("clothoLocalStorage",["$window","PubSub","Debug","$document",function(a,b,c,d){function e(){function d(a,c){b.trigger("update",a),b.trigger("update:"+a,angular.copy(c))}var e=new c("clothoLocalStorage","#88cc88"),f=a.localStorage,g=JSON,h="clotho_",i=function(){try{if(f.setItem("test","7"),"7"===f.getItem("test"))return f.removeItem("test"),!0}catch(a){}return!1};e.log("localStorage support? "+i());var j=function(){for(var a=!1,b=0,c=f.length;c>b;++b){var d=f.key(b);angular.isString(d)&&d.substring(0,h.length)==h&&(f.removeItem(d),--b,--c,a=!0)}return a},k=function(a,b){if(angular.isEmpty(a))return null;var c=f.getItem(h+a);return angular.isEmpty(c)?angular.isDefined(b)?b:null:angular.isObject(c)?c:g.parse(c)},l=function(a){if(angular.isEmpty(a))return!1;var b=f.getItem(h+a);return angular.isDefined(b)&&!angular.isEmpty(b)},m=function(a){return angular.isEmpty(a)?!1:(f.removeItem(h+a),!0)},n=function(a,b){if(angular.isEmpty(a))return!1;if(angular.isEmpty(b))return!1;b=angular.isString(b)?b:g.stringify(b);try{f.setItem(h+a,b)}catch(c){e.error("couldnt save "+a,c)}return!0},o=function(b){if(b||(b=a.event),!(b.key.indexOf(h)<0)){var c=b.key.replace(h,"")||"",f=k(c);e.log("handle_storage_event for "+c),d(c,f)}};return a.addEventListener?a.addEventListener("storage",o,!1):a.attachEvent("onstorage",o),{getPrefix:function(){return h},isSupported:i,clear:j,getItem:k,hasItem:l,removeItem:m,setItem:n}}return a.$clotho.$localStorage?a.$clotho.$localStorage:a.$clotho.$localStorage=e()}]),angular.module("clotho.core").service("Collector",["$window","clothoLocalStorage","PubSub","Debug",function(a,b,c,d){function e(){function a(a,b){c.trigger("update",a),c.trigger("update:"+a,angular.copy(b))}var e=new d("Collector","#55bb55"),f=function(a,c){b.setItem(a,c)},g=function(b,c,d){d||!angular.equals(i(b),c)?(e.log(b+" (saving)",c),f(b,c),a(b,c)):e.log(b+" (model unchanged)")},h=function(a){return angular.copy(i(a))},i=function(a){return b.hasItem(a)?b.getItem(a):!1},j=function(a){return b.hasItem(a)?(b.removeItem(a),!0):!1},k=function(){b.clear(),c.trigger("collector_reset")};return{hasItem:b.hasItem,silentAddModel:f,storeModel:g,retrieveModel:h,retrieveRef:i,removeModel:j,clearStorage:k}}return a.$clotho.$collector?a.$clotho.$collector:a.$clotho.$collector=e()}]),angular.module("clotho.core").service("PubSub",["$window","$rootScope","$filter","Debug",function(a,b,c,d){function e(){function a(a){return a?"all"==a?Object.keys(k):a.split(l):[]}function e(a){return k.hasOwnProperty(a)}function f(a){return e(a)&&k[a].length>0}function g(a,b,c,d){return angular.isFunction(a)?{callback:a,ref:b||null,priority:parseInt(c)||100,once:1==d}:null}function h(a,b){if(b&&!angular.isEmpty(b)){k[a]||(k[a]=[]);var c=b.ref;return angular.isScope(c)&&c.$on("$destroy",function(){t(j(c))}),k[a].push(b),angular.once(function(){i(a,b)})}return angular.noop}function i(a,b){var c=angular.remove(k[a],function(a){return angular.equals(a,b)});return c.length>0}function j(a){return angular.isScope(a)?a.$id:a}var k={},l=/\s+/,m=new d("PubSub","#55cccc"),n=function(){m.log("LISTENERS:"),angular.forEach(k,function(a,b){m.log(b),m.table(k[b])})},o=function(d,e){e=angular.isUndefined(e)||null===e?null:[e],angular.forEach(a(d),function(a){m.log("Publish on "+a,e),f(a)&&angular.forEach(c("orderBy")(k[a],"priority"),function(c,d){b.$safeApply(function(){c.callback.apply(c.ref,e)}),1==c.once&&k[a].splice(d,1)})})},p=function(c){angular.forEach(a(c),function(a){m.log("Reject on "+a),angular.forEach(k[a],function(c,d){b.$safeApply(function(){c.callback.apply(c.ref,null)}),1==c.once&&k[a].splice(d,1)})})},q=function(b,c,d,e,f){d=d||null,f=1==f,f&&(c=angular.once(c));var i=[];return angular.forEach(a(b),function(a){i.push(h(a,g(c,d,e,f)))}),function(){angular.forEach(i,function(a){a()})}},r=function(a,b,c,d){q(a,b,c,d,!0)},s=function(a,b){var c=angular.remove(k[a],function(a){return angular.equals(a.callback,b)});return c.length>0},t=function(a){angular.forEach(k,function(b){angular.remove(b,function(b){return angular.equals(j(a),j(b.ref))})})},u=function(b){angular.forEach(a(b),function(a,b){k[b].length=0})};return{logListeners:n,trigger:o,on:q,once:r,off:s,destroy:t,reject:p,clear:u}}return a.$clotho.$pubsub?a.$clotho.$pubsub:a.$clotho.$pubsub=e()}]),angular.module("clotho.core").service("Socket",["$window","$q","$log","PubSub","ClientAPI","Debug",function(a,b,c,d,e,f){function g(){function b(a){return h?(a=angular.isObject(a)?JSON.stringify(a):a,j.log("sending data: "+a),void g.send(a)):(j.log("(not ready) queueing request: ",a),void i.push(a))}function c(){angular.forEach(i,function(a){b(a)}),i=[]}var g,h,i=[],j=new f("Socket","#5555bb");return g=a.$clotho.socket?a.$clotho.socket:a.$clotho.socket=new WebSocket("wss://"+window.location.host+window.location.pathname+"websocket"),1==g.readyState?(j.log("already exists, sending items in socket Queue..."),h=!0,c()):g.onopen=function(){j.log("opened, sending queued items..."),h=!0,c()},g.onerror=function(a){j.error("socket error",a)},g.onclose=function(){h=!1,e.say({"class":"error",text:"Socket Connection Closed",from:"client"})},g.onmessage=function(a){a=JSON.parse(a.data),j.log("received",a);var b=a.channel,c=a.requestId,f=angular.isUndefined(a.data),g=a.data;if(angular.isFunction(e[b]))e[b](g);else{var h=f?"reject":"trigger";c?d[h](b+":"+c,g):d[h](b,g)}},{state:function(){return g.readyState},emit:function(a,c,d){var e={channel:a,data:c};b(e),"function"==typeof d&&d(e)},send:b}}return a.$clotho.$socket?a.$clotho.$socket:a.$clotho.$socket=g()}]);
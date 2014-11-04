angular.module("clotho.interface",["ui.bootstrap","ui.bootstrap-decorate","ui.codemirror","cfp.hotkeys","angularFileUpload","chieffancypants.loadingBar","clotho.urlShorten"]),angular.module("clotho.interface").service("interfaceConfig",["$injector",function(a){return{loginVisible:a.has("clothoEditorDirective")}}]),angular.module("ui.codemirror",[]).constant("uiCodemirrorConfig",{}).directive("uiCodemirror",["uiCodemirrorConfig",function(a){return{restrict:"EA",require:"?ngModel",priority:1,compile:function(){if(angular.isUndefined(window.CodeMirror))throw new Error("ui-codemirror need CodeMirror to work... (o rly?)");return function(b,c,d,e){var f,g,h,i;if(i=c.text(),f=a.codemirror||{},g=angular.extend({value:i},f,b.$eval(d.uiCodemirror),b.$eval(d.uiCodemirrorOpts)),"TEXTAREA"===c[0].tagName?h=window.CodeMirror.fromTextArea(c[0],g):(c.html(""),h=new window.CodeMirror(function(a){c.append(a)},g)),d.uiCodemirror||d.uiCodemirrorOpts){var j=Object.keys(window.CodeMirror.defaults);b.$watch(d.uiCodemirror||d.uiCodemirrorOpts,function(a,b){angular.isObject(a)&&j.forEach(function(c){if(a.hasOwnProperty(c)){if(b&&a[c]===b[c])return;h.setOption(c,a[c])}})},!0)}e&&(e.$formatters.push(function(a){if(angular.isUndefined(a)||null===a)return"";if(angular.isObject(a)||angular.isArray(a))throw new Error("ui-codemirror cannot use an object or an array as a model");return a}),e.$render=function(){var a=e.$viewValue||"";h.setValue(a)},h.on("change",function(a){var c=a.getValue();c!==e.$viewValue&&b.$apply(function(){e.$setViewValue(c)})})),d.uiRefresh&&b.$watch(d.uiRefresh,function(a,b){a!==b&&h.refresh()}),b.$on("CodeMirror",function(a,b){if(!angular.isFunction(b))throw new Error("the CodeMirror event requires a callback function");b(h)}),angular.isFunction(g.onLoad)&&g.onLoad(h)}}}}]),function(){var a=angular.module("angularFileUpload",[]);a.service("$upload",["$http","$timeout",function(a,b){function c(c){c.method=c.method||"POST",c.headers=c.headers||{},c.transformRequest=c.transformRequest||function(b,c){return window.ArrayBuffer&&b instanceof window.ArrayBuffer?b:a.defaults.transformRequest[0](b,c)},window.XMLHttpRequest.__isShim&&(c.headers.__setXHR_=function(){return function(a){c.__XHR=a,c.xhrFn&&c.xhrFn(a),a.upload.addEventListener("progress",function(a){c.progress&&b(function(){c.progress&&c.progress(a)})},!1),a.upload.addEventListener("load",function(a){a.lengthComputable&&b(function(){c.progress&&c.progress(a)})},!1)}});var d=a(c);return d.progress=function(a){return c.progress=a,d},d.abort=function(){return c.__XHR&&b(function(){c.__XHR.abort()}),d},d.xhr=function(a){return c.xhrFn=a,d},d.then=function(a,b){return function(d,e,f){c.progress=f||c.progress;var g=b.apply(a,[d,e,f]);return g.abort=a.abort,g.progress=a.progress,g.xhr=a.xhr,g}}(d,d.then),d}this.upload=function(b){b.headers=b.headers||{},b.headers["Content-Type"]=void 0,b.transformRequest=b.transformRequest||a.defaults.transformRequest;var d=new FormData,e=b.transformRequest,f=b.data;return b.transformRequest=function(a,c){if(f)if(b.formDataAppender)for(var d in f){var g=f[d];b.formDataAppender(a,d,g)}else for(var d in f){var g=f[d];if("function"==typeof e)g=e(g,c);else for(var h=0;h<e.length;h++){var i=e[h];"function"==typeof i&&(g=i(g,c))}a.append(d,g)}if(null!=b.file){var j=b.fileFormDataName||"file";if("[object Array]"===Object.prototype.toString.call(b.file))for(var k="[object String]"===Object.prototype.toString.call(j),h=0;h<b.file.length;h++)a.append(k?j+h:j[h],b.file[h],b.file[h].name);else a.append(j,b.file,b.file.name)}return a},b.data=d,c(b)},this.http=function(a){return c(a)}}]),a.directive("ngFileSelect",["$parse","$timeout",function(a,b){return function(c,d,e){var f=a(e.ngFileSelect);d.bind("change",function(a){var d,e,g=[];if(d=a.target.files,null!=d)for(e=0;e<d.length;e++)g.push(d.item(e));b(function(){f(c,{$files:g,$event:a})})}),d.bind("click",function(){this.value=null})}}]),a.directive("ngFileDropAvailable",["$parse","$timeout",function(a,b){return function(c,d,e){if("draggable"in document.createElement("span")){var f=a(e.ngFileDropAvailable);b(function(){f(c)})}}}]),a.directive("ngFileDrop",["$parse","$timeout",function(a,b){return function(c,d,e){if("draggable"in document.createElement("span")){var f=null,g=a(e.ngFileDrop);d[0].addEventListener("dragover",function(a){b.cancel(f),a.stopPropagation(),a.preventDefault(),d.addClass(e.ngFileDragOverClass||"dragover")},!1),d[0].addEventListener("dragleave",function(){f=b(function(){d.removeClass(e.ngFileDragOverClass||"dragover")})},!1),d[0].addEventListener("drop",function(a){a.stopPropagation(),a.preventDefault(),d.removeClass(e.ngFileDragOverClass||"dragover");var f,h=[],i=a.dataTransfer.files;if(null!=i)for(f=0;f<i.length;f++)h.push(i.item(f));b(function(){g(c,{$files:h,$event:a})})},!1)}}}])}(),function(){"use strict";angular.module("cfp.hotkeys",[]).provider("hotkeys",function(){this.includeCheatSheet=!0,this.templateTitle="Keyboard Shortcuts:",this.template='<div class="cfp-hotkeys-container fade" ng-class="{in: helpVisible}" style="display: none;"><div class="cfp-hotkeys"><h4 class="cfp-hotkeys-title">{{ title }}</h4><table><tbody><tr ng-repeat="hotkey in hotkeys | filter:{ description: \'!$$undefined$$\' }"><td class="cfp-hotkeys-keys"><span ng-repeat="key in hotkey.format() track by $index" class="cfp-hotkeys-key">{{ key }}</span></td><td class="cfp-hotkeys-text">{{ hotkey.description }}</td></tr></tbody></table><div class="cfp-hotkeys-close" ng-click="toggleCheatSheet()">×</div></div></div>',this.cheatSheetHotkey="?",this.cheatSheetDescription="Show / hide this help menu",this.$get=["$rootElement","$rootScope","$compile","$window","$document",function(a,b,c,d,e){function f(a){var b={command:"⌘",shift:"⇧",left:"←",right:"→",up:"↑",down:"↓","return":"↩",backspace:"⌫"};a=a.split("+");for(var c=0;c<a.length;c++)"mod"===a[c]&&(a[c]=d.navigator&&d.navigator.platform.indexOf("Mac")>=0?"command":"ctrl"),a[c]=b[a[c]]||a[c];return a.join(" + ")}function g(a,b,c,d,e,f){this.combo=a instanceof Array?a:[a],this.description=b,this.callback=c,this.action=d,this.allowIn=e,this.persistent=f}function h(){for(var a=o.hotkeys.length;a--;){var b=o.hotkeys[a];b&&!b.persistent&&k(b)}}function i(){o.helpVisible=!o.helpVisible,o.helpVisible?(t=l("esc"),k("esc"),j("esc",t.description,i)):(k("esc"),t!==!1&&j(t))}function j(a,b,c,d,e,f){var h,i=["INPUT","SELECT","TEXTAREA"],j=Object.prototype.toString.call(a);if("[object Object]"===j&&(b=a.description,c=a.callback,d=a.action,f=a.persistent,e=a.allowIn,a=a.combo),b instanceof Function?(d=c,c=b,b="$$undefined$$"):angular.isUndefined(b)&&(b="$$undefined$$"),void 0===f&&(f=!0),"function"==typeof c){h=c,e instanceof Array||(e=[]);for(var k,l=0;l<e.length;l++)e[l]=e[l].toUpperCase(),k=i.indexOf(e[l]),-1!==k&&i.splice(k,1);c=function(a){var b=!0,c=a.target||a.srcElement,d=c.nodeName.toUpperCase();if((" "+c.className+" ").indexOf(" mousetrap ")>-1)b=!0;else for(var e=0;e<i.length;e++)if(i[e]===d){b=!1;break}b&&n(h.apply(this,arguments))}}"string"==typeof d?Mousetrap.bind(a,n(c),d):Mousetrap.bind(a,n(c));var m=new g(a,b,c,d,e,f);return o.hotkeys.push(m),m}function k(a){var b=a instanceof g?a.combo:a;if(Mousetrap.unbind(b),b instanceof Array){for(var c=!0,d=0;d<b.length;d++)c=k(b[d])&&c;return c}var e=o.hotkeys.indexOf(l(b));return e>-1?(o.hotkeys[e].combo.length>1?o.hotkeys[e].combo.splice(o.hotkeys[e].combo.indexOf(b),1):o.hotkeys.splice(e,1),!0):!1}function l(a){for(var b,c=0;c<o.hotkeys.length;c++)if(b=o.hotkeys[c],b.combo.indexOf(a)>-1)return b;return!1}function m(a){return p[a.$id]=[],a.$on("$destroy",function(){for(var b=p[a.$id].length;b--;)k(p[a.$id][b]),delete p[a.$id][b]}),{add:function(b){var c;return c=arguments.length>1?j.apply(this,arguments):j(b),p[a.$id].push(c),this}}}function n(a){return function(c,d){if(a instanceof Array){var e=a[0],f=a[1];a=function(){f.scope.$eval(e)}}b.$apply(function(){a(c,l(d))})}}Mousetrap.stopCallback=function(a,b){return(" "+b.className+" ").indexOf(" mousetrap ")>-1?!1:b.contentEditable&&"true"==b.contentEditable},g.prototype.format=function(){for(var a=this.combo[0],b=a.split(/[\s]/),c=0;c<b.length;c++)b[c]=f(b[c]);return b};var o=b.$new();o.hotkeys=[],o.helpVisible=!1,o.title=this.templateTitle,o.toggleCheatSheet=i;var p=[];if(b.$on("$routeChangeSuccess",function(a,b){h(),b&&b.hotkeys&&angular.forEach(b.hotkeys,function(a){var c=a[2];("string"==typeof c||c instanceof String)&&(a[2]=[c,b]),a[5]=!1,j.apply(this,a)})}),this.includeCheatSheet){var q=e[0],r=a[0],s=angular.element(this.template);j(this.cheatSheetHotkey,this.cheatSheetDescription,i),(r===q||r===q.documentElement)&&(r=q.body),angular.element(r).append(c(s)(o))}var t=!1,u={add:j,del:k,get:l,bindTo:m,template:this.template,toggleCheatSheet:i,includeCheatSheet:this.includeCheatSheet,cheatSheetHotkey:this.cheatSheetHotkey,cheatSheetDescription:this.cheatSheetDescription,purgeHotkeys:h,templateTitle:this.templateTitle};return u}]}).directive("hotkey",["hotkeys",function(a){return{restrict:"A",link:function(b,c,d){var e,f;angular.forEach(b.$eval(d.hotkey),function(b,c){f="string"==typeof d.hotkeyAllowIn?d.hotkeyAllowIn.split(/[\s,]+/):[],e=c,a.add({combo:c,description:d.hotkeyDescription,callback:b,action:d.hotkeyAction,allowIn:f})}),c.bind("$destroy",function(){a.del(e)})}}}]).run(["hotkeys",function(){}])}(),function(a,b){function c(a,b,c){return a.addEventListener?void a.addEventListener(b,c,!1):void a.attachEvent("on"+b,c)}function d(a){if("keypress"==a.type){var b=String.fromCharCode(a.which);return a.shiftKey||(b=b.toLowerCase()),b}return y[a.which]?y[a.which]:z[a.which]?z[a.which]:String.fromCharCode(a.which).toLowerCase()}function e(a,b){return a.sort().join(",")===b.sort().join(",")}function f(a){a=a||{};var b,c=!1;for(b in E)a[b]?c=!0:E[b]=0;c||(H=!1)}function g(a,b,c,d,f,g){var h,i,j=[],k=c.type;if(!C[a])return[];for("keyup"==k&&n(a)&&(b=[a]),h=0;h<C[a].length;++h)if(i=C[a][h],(d||!i.seq||E[i.seq]==i.level)&&k==i.action&&("keypress"==k&&!c.metaKey&&!c.ctrlKey||e(b,i.modifiers))){var l=!d&&i.combo==f,m=d&&i.seq==d&&i.level==g;(l||m)&&C[a].splice(h,1),j.push(i)}return j}function h(a){var b=[];return a.shiftKey&&b.push("shift"),a.altKey&&b.push("alt"),a.ctrlKey&&b.push("ctrl"),a.metaKey&&b.push("meta"),b}function i(a){return a.preventDefault?void a.preventDefault():void(a.returnValue=!1)}function j(a){return a.stopPropagation?void a.stopPropagation():void(a.cancelBubble=!0)}function k(a,b,c,d){J.stopCallback(b,b.target||b.srcElement,c,d)||a(b,c)===!1&&(i(b),j(b))}function l(a,b,c){var d,e=g(a,b,c),h={},i=0,j=!1;for(d=0;d<e.length;++d)e[d].seq&&(i=Math.max(i,e[d].level));for(d=0;d<e.length;++d)if(e[d].seq){if(e[d].level!=i)continue;j=!0,h[e[d].seq]=1,k(e[d].callback,c,e[d].combo,e[d].seq)}else j||k(e[d].callback,c,e[d].combo);var l="keypress"==c.type&&G;c.type!=H||n(a)||l||f(h),G=j&&"keydown"==c.type}function m(a){"number"!=typeof a.which&&(a.which=a.keyCode);var b=d(a);if(b)return"keyup"==a.type&&F===b?void(F=!1):void J.handleKey(b,h(a),a)}function n(a){return"shift"==a||"ctrl"==a||"alt"==a||"meta"==a}function o(){clearTimeout(x),x=setTimeout(f,1e3)}function p(){if(!w){w={};for(var a in y)a>95&&112>a||y.hasOwnProperty(a)&&(w[y[a]]=a)}return w}function q(a,b,c){return c||(c=p()[a]?"keydown":"keypress"),"keypress"==c&&b.length&&(c="keydown"),c}function r(a,b,c,e){function g(b){return function(){H=b,++E[a],o()}}function h(b){k(c,b,a),"keyup"!==e&&(F=d(b)),setTimeout(f,10)}E[a]=0;for(var i=0;i<b.length;++i){var j=i+1===b.length,l=j?h:g(e||t(b[i+1]).action);u(b[i],l,e,a,i)}}function s(a){return"+"===a?["+"]:a.split("+")}function t(a,b){var c,d,e,f=[];for(c=s(a),e=0;e<c.length;++e)d=c[e],B[d]&&(d=B[d]),b&&"keypress"!=b&&A[d]&&(d=A[d],f.push("shift")),n(d)&&f.push(d);return b=q(d,f,b),{key:d,modifiers:f,action:b}}function u(a,b,c,d,e){D[a+":"+c]=b,a=a.replace(/\s+/g," ");var f,h=a.split(" ");return h.length>1?void r(a,h,b,c):(f=t(a,c),C[f.key]=C[f.key]||[],g(f.key,f.modifiers,{type:f.action},d,a,e),void C[f.key][d?"unshift":"push"]({callback:b,modifiers:f.modifiers,action:f.action,seq:d,level:e,combo:a}))}function v(a,b,c){for(var d=0;d<a.length;++d)u(a[d],b,c)}for(var w,x,y={8:"backspace",9:"tab",13:"enter",16:"shift",17:"ctrl",18:"alt",20:"capslock",27:"esc",32:"space",33:"pageup",34:"pagedown",35:"end",36:"home",37:"left",38:"up",39:"right",40:"down",45:"ins",46:"del",91:"meta",93:"meta",224:"meta"},z={106:"*",107:"+",109:"-",110:".",111:"/",186:";",187:"=",188:",",189:"-",190:".",191:"/",192:"`",219:"[",220:"\\",221:"]",222:"'"},A={"~":"`","!":"1","@":"2","#":"3",$:"4","%":"5","^":"6","&":"7","*":"8","(":"9",")":"0",_:"-","+":"=",":":";",'"':"'","<":",",">":".","?":"/","|":"\\"},B={option:"alt",command:"meta","return":"enter",escape:"esc",mod:/Mac|iPod|iPhone|iPad/.test(navigator.platform)?"meta":"ctrl"},C={},D={},E={},F=!1,G=!1,H=!1,I=1;20>I;++I)y[111+I]="f"+I;for(I=0;9>=I;++I)y[I+96]=I;c(b,"keypress",m),c(b,"keydown",m),c(b,"keyup",m);var J={bind:function(a,b,c){return a=a instanceof Array?a:[a],v(a,b,c),this},unbind:function(a,b){return J.bind(a,function(){},b)},trigger:function(a,b){return D[a+":"+b]&&D[a+":"+b]({},a),this},reset:function(){return C={},D={},this},stopCallback:function(a,b){return(" "+b.className+" ").indexOf(" mousetrap ")>-1?!1:"INPUT"==b.tagName||"SELECT"==b.tagName||"TEXTAREA"==b.tagName||b.isContentEditable},handleKey:l};a.Mousetrap=J,"function"==typeof define&&define.amd&&define(J)}(window,document),function(){"use strict";angular.module("angular-loading-bar",["cfp.loadingBarInterceptor"]),angular.module("chieffancypants.loadingBar",["cfp.loadingBarInterceptor"]),angular.module("cfp.loadingBarInterceptor",["cfp.loadingBar"]).config(["$httpProvider",function(a){var b=["$q","$cacheFactory","$timeout","$rootScope","cfpLoadingBar",function(b,c,d,e,f){function g(){d.cancel(i),f.complete(),k=0,j=0}function h(b){var d,e=a.defaults;if("GET"!==b.method||b.cache===!1)return b.cached=!1,!1;d=b.cache===!0&&void 0===e.cache?c.get("$http"):void 0!==e.cache?e.cache:b.cache;var f=void 0!==d?void 0!==d.get(b.url):!1;return void 0!==b.cached&&f!==b.cached?b.cached:(b.cached=f,f)}var i,j=0,k=0,l=f.latencyThreshold;return{request:function(a){return a.ignoreLoadingBar||h(a)||(e.$broadcast("cfpLoadingBar:loading",{url:a.url}),0===j&&(i=d(function(){f.start()},l)),j++,f.set(k/j)),a},response:function(a){return a.config.ignoreLoadingBar||h(a.config)||(k++,e.$broadcast("cfpLoadingBar:loaded",{url:a.config.url}),k>=j?g():f.set(k/j)),a},responseError:function(a){return a.config.ignoreLoadingBar||h(a.config)||(k++,e.$broadcast("cfpLoadingBar:loaded",{url:a.config.url}),k>=j?g():f.set(k/j)),b.reject(a)}}}];a.interceptors.push(b)}]),angular.module("cfp.loadingBar",[]).provider("cfpLoadingBar",function(){this.includeSpinner=!0,this.includeBar=!0,this.latencyThreshold=100,this.startSize=.02,this.parentSelector="body",this.$get=["$document","$timeout","$animate","$rootScope",function(a,b,c,d){function e(){var e=a.find(l);b.cancel(k),p||(d.$broadcast("cfpLoadingBar:started"),p=!0,s&&c.enter(m,e),r&&c.enter(o,e),f(t))}function f(a){if(p){var c=100*a+"%";n.css("width",c),q=a,b.cancel(j),j=b(function(){g()},250)}}function g(){if(!(h()>=1)){var a=0,b=h();a=b>=0&&.25>b?(3*Math.random()+3)/100:b>=.25&&.65>b?3*Math.random()/100:b>=.65&&.9>b?2*Math.random()/100:b>=.9&&.99>b?.005:0;var c=h()+a;f(c)}}function h(){return q}function i(){d.$broadcast("cfpLoadingBar:completed"),f(1),k=b(function(){c.leave(m,function(){q=0,p=!1}),c.leave(o)},500)}var j,k,l=this.parentSelector,m=angular.element('<div id="loading-bar"><div class="bar"><div class="peg"></div></div></div>'),n=m.find("div").eq(0),o=angular.element('<div id="loading-bar-spinner"><div class="spinner-icon"></div></div>'),p=!1,q=0,r=this.includeSpinner,s=this.includeBar,t=this.startSize;return{start:e,set:f,status:h,inc:g,complete:i,includeSpinner:this.includeSpinner,latencyThreshold:this.latencyThreshold,parentSelector:this.parentSelector,startSize:this.startSize}}]})}(),angular.module("clotho.interface").directive("formField",["$parse",function(a){var b='<ng-form class="form-group" ng-transclude></ng-form>';return{restrict:"E",template:b,replace:!0,transclude:!0,require:["^form"],controller:["$scope","$element","$attrs",function(){}],link:function(b,c,d,e){var f=e[0],g=d.name,h=d.help,i=a(d.horizontal)(b),j=(d.removableField,!1),k=c.children();if(1!==k.length)throw"You must include the form input element as the single child of this directive, instead got "+c.html()+"generating:";var l=k.attr("id"),m=(k.prop("tagName"),k.attr("type")),n=/checkbox|radio|button/gi,o=/file|checkbox|radio/gi;l||(l="input"+Math.floor(1e7*Math.random()).toString(),k.attr("id",l)),o.test(m)||angular.isDefined(d.noStyling)||k.addClass("form-control");var p=angular.element('<label class="control-label" for="'+l+'">'+("radio"===m?k.val():g)+"</label>");if(n.test(m)){var q=angular.element('<div class="'+m+'"></div>');q.append(k),k=q}if(i){p.addClass("col-sm-3");var q=angular.element('<div class="col-sm-9"></div>');q.append(k)}c.append(q),g&&!j&&c.prepend(p),h&&c.append('<p class="help-block">'+h+"</p>"),b.$watch(function(){return f.$invalid},function(a){c.toggleClass("has-error",a)})}}}]),angular.module("ui.bootstrap-decorate",["ui.bootstrap"]).config(["$provide",function(a){a.decorator("$modal",["$delegate",function(a){return a.messageBox=function(b,c,d){return a.open({templateUrl:"views/_interface/ui-custom/dialogMessagebox.html",controller:"MessageBoxController",resolve:{model:function(){return{title:b,message:c,buttons:d}}}})},a.share=function(b){return a.open({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,templateUrl:"views/_interface/ui-custom/dialogShare.html",controller:"DialogShareController",resolve:{model:function(){return{url:b}}}})},a.video=function(b,c){return angular.extend(c,{width:"560"}),a.open({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,template:'<div youtube="{{ videoId }}" params="videoParams""></div>',controller:"VideoDialogController",resolve:{model:function(){return{videoId:b,videoParams:c}}}})},a.clothoEdit=function(b){return a.open({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!1,templateUrl:"views/_interface/ui-custom/clothoEditModal.html",controller:"DialogEditController",resolve:{model:function(){return{uuid:b}}}})},a}])}]),angular.module("ui.bootstrap-decorate").controller("MessageBoxController",["$scope","$modalInstance","model",function(a,b,c){a.title=c.title,a.message=c.message,a.buttons=c.buttons,a.close=function(a){b.close(a)}}]),angular.module("ui.bootstrap-decorate").controller("DialogShareController",["$scope","$modalInstance","model","$location","$window",function(a,b){a.close=function(a){b.close(a)}}]),angular.module("ui.bootstrap-decorate").controller("VideoDialogController",["$scope","$modalInstance","model",function(a,b,c){a.videoId=c.videoId,a.videoParams=c.videoParams,a.close=function(a){b.close(a)}}]),angular.module("ui.bootstrap-decorate").controller("DialogEditController",["$scope","$modalInstance","model",function(a,b,c){a.uuid=c.uuid}]),angular.module("clotho.interface").directive("popoverHtmlUnsafePopup",["$templateCache",function(a){return a.put("template/popover/popover-html-unsafe-popup.html",'<div class="popover {{placement}}" ng-class="{ in: isOpen(), fade: animation() }">\n  <div class="arrow"></div>\n\n  <div class="popover-inner">\n      <h3 class="popover-title" ng-bind="title" ng-show="title"></h3>\n      <div class="popover-content" bind-html-unsafe="content"></div>\n  </div>\n</div>\n'),{restrict:"EA",replace:!0,scope:{title:"@",content:"@",placement:"@",animation:"&",isOpen:"&"},templateUrl:"template/popover/popover-html-unsafe-popup.html"}}]).directive("popoverHtmlUnsafe",["$compile","$timeout","$parse","$window","$tooltip",function(a,b,c,d,e){return e("popoverHtmlUnsafe","popover","click")}]),angular.module("clotho.interface").directive("clothoSharable",["$injector","Clotho","ClothoSchemas",function(a,b,c){var d=a.has("clothoEditorDirective");return{restrict:"A",scope:{sharable:"=clothoSharable"},templateUrl:"views/_interface/clothoSharable.html",link:function(a){a.$watch("sharable",function(b){a.type=c.dirtyDetermineType(b),a.iconClass=c.determineSharableIcon(a.type),a.labelClass="label-"+c.typeToColorClass(a.type),a.schemaName=c.mapSchemaIdToName(b.schema)}),a.editorPresent=d,a.edit=b.edit}}}]),angular.module("clotho.clothoDirectives").directive("restrictionEnzymePopup",["$clothoPopup",function(a){return a("restrictionEnzymePopup",{placement:"top",trigger:"mouseenter"})}]).directive("restrictionEnzymePopupInner",function(){return{restrict:"A",replace:!0,scope:{enzyme:"=?sharableModel",placement:"@popupPlacement",reposition:"&"},templateUrl:"views/_interface/dna/restrictionEnzymePopup.html"}}),angular.module("clotho.clothoDirectives").directive("sequencePopup",["$clothoPopup",function(a){return a("sequencePopup",{placement:"top",trigger:"mouseenter"})}]).directive("sequencePopupInner",function(){return{restrict:"A",replace:!0,scope:{inputModel:"=?sharableModel",title:"=popupTitle",placement:"@popupPlacement",reposition:"&"},templateUrl:"views/_interface/dna/sequencePopup.html",link:function(a){function b(b,c){return">"+(c||b.name||a.title)+(b.description?" "+b.description:"")+"\n"+b.sequence}a.$watch("inputModel",function(c){if(a.model=angular.isArray(c)?c:[c],c){var d=angular.map(a.model,function(c,d){var e=c.name||(a.title?a.title+" - ":"")+"Fragment "+(d+1);return b(c,e)});a.generatedFASTA=d.join("\n\n")}})}}}),angular.module("clotho.interface").directive("textfileDownload",["$window","$document",function(a,b){const c="text/plain";return{restrict:"E",scope:{content:"=",name:"@?"},link:function(d,e){function f(){var f=d.content,g=d.name;if(!angular.isUndefined(f)){var h=new a.Blob([f],{type:c}),i=b[0].createElement("a");i.download=g,i.alt="Download as File",i.href=a.URL.createObjectURL(h),i.innerHTML='<span class="glyphicon glyphicon-download"></span>',e.replaceWith(i)}}angular.isUndefined(a.Blob)||angular.isUndefined(a.URL.createObjectURL)||(d.$watch("content",f),d.$watch("name",f))}}}]),angular.module("clotho.interface").directive("closeable",["$compile",function(a){return{restrict:"A",priority:1100,link:function(b,c){b.removeElement=function(){c.remove()},c.prepend(a('<a class="close" style="position: absolute; top: 12px; right: 15px;" ng-click="removeElement()">&times;</a>')(b))}}}]),angular.module("clotho.interface").directive("hintButton",function(){return{replace:!0,template:'<button class="btn" popover="{{ hint }}" popover-trigger="mouseenter" popover-placement="left"><i class="glyphicon glyphicon-info-sign"></i> Hint</button>',link:function(a,b,c){a.hint=c.hintButton}}}),angular.module("clotho.interface").directive("uiEvent",["$parse",function(a){return function(b,c,d){var e=b.$eval(d.uiEvent);angular.forEach(e,function(d,e){var f=a(d);c.bind(e,function(a){var c=Array.prototype.slice.call(arguments);c=c.splice(1),f(b,{$event:a,$params:c}),b.$$phase||b.$apply()})})}}]),angular.module("clotho.interface").directive("restrictInput",function(){return{restrict:"A",require:"ngModel",link:function(a,b,c,d){var e=function(b){var e=a.$eval(c.restrictInput),f=b.replace(e,"");return f!=b&&(d.$setViewValue(f),d.$render()),f};d.$parsers.unshift(e)}}}),angular.module("clotho.interface").directive("clickOutside",["$document","$parse",function(a,b){return function(c,d,e){var f,g=b(e.clickOutside);c.$watch(function(){return b(e.clickOutsideActive)(c)},function(a){f=!!a});var h=function(a){f&&(d[0].contains(a.target)||(a.preventDefault(),a.stopPropagation(),c.$apply(g(c,{$event:a}))))};a.bind("click",h),c.$on("$destroy",function(){a.unbind("click",h)})}}]),angular.module("clotho.interface").directive("simpleSticky",["$window",function(a){function b(){return angular.isDefined(a.pageYOffset)?a.pageYOffset:a.document[0].documentElement.scrollTop}function c(a,b){f.push(a),b.$on("$destroy",d)}function d(a){_.remove(f,a)}var e=angular.element(a),f=[],g=_.throttle(function(){var a=b();angular.forEach(f,function(b){b.apply(null,[a])})},50);return e.on("scroll resize",g),{restrict:"A",link:function(a,d,e){function f(a){var b=a+i>l;b!==h&&(h=b,g(b))}function g(a){d.css(a?{top:i+"px",width:"inherit",position:"fixed"}:{top:k,position:j})}var h,i=parseInt(e.simpleSticky,10)||70,j=d.css("position"),k=d.css("top"),l=d[0].getBoundingClientRect().top+b();c(function(a){f(a)},a),f()}}}]),angular.module("clotho.interface").filter("breakLines",function(){return function(a,b,c){return(a.match(new RegExp(".{1,"+b+"}","gi"))||[]).join(c||"\n")}}),angular.module("clotho.interface").filter("categorize",["$parse",function(a){return function(b,c){if(c===!1)return b;if((c||angular.isUndefined(c))&&angular.isArray(b)){var d={},e=angular.isString(c)?a(c):function(a){return a},f=function(a){return angular.isObject(a)?e(a):a};angular.forEach(b,function(a){var b=f(a);d[b]||(d[b]=[]),d[b].push(a)}),b=d}return b}}]),angular.module("clotho.interface").filter("highlight",function(){return function(a,b,c,d){return b||angular.isNumber(b)?(c=angular.isDefined(c)?c:"ui-match",a=a.toString(),b=b.toString(),d?a.split(b).join('<span class="'+c+'">'+b+"</span>"):a.replace(new RegExp(b,"gi"),'<span class="'+c+'">$&</span>')):a}}),angular.module("clotho.interface").filter("joinArray",function(){return function(a,b){if(angular.isObject(a)){if(angular.isArray(a))return a.join(b);var c="";return angular.forEach(a,function(a,d){c+=d+b}),c.substring(0,c.length-b.length)}return a}}),angular.module("clotho.interface").filter("plainText",function(){return function(a,b){return b?a.replace(/(<([^>]+)>)/gi,""):a}}),angular.module("clotho.interface").filter("shuffle",function(){return function(a,b){if(b!==!1&&angular.isArray(items)){for(var c,d,e=a.length;e;)d=Math.floor(Math.random()*e--),c=a[e],a[e]=a[d],a[d]=c;return a}return items}}),angular.module("clotho.interface").filter("stringEnds",function(){return function(a,b){b=b||20;var c=function(a){return a.substr(0,b)+"..."+a.substring(a.length-b)};if(!angular.isString(a)){if(angular.isArray(a)){var d=angular.isString(a[0])?c(a[0]):"";return"<Array>"+(d?" First Entry:\n"+d:"")}if(angular.isObject(a)){var e=angular.isString(a.key)?c(a.key):"";return"<Object>"+(e?"Object.key:\n"+e:"")}return a}return a.length<=2*b?a:c(a)}}),angular.module("clotho.interface").filter("unique",["$parse",function(a){return function(b,c){if(c===!1)return b;if((c||angular.isUndefined(c))&&angular.isArray(b)){var d=[],e=angular.isString(c)?a(c):function(a){return a},f=function(a){return angular.isObject(a)?e(a):a};angular.forEach(b,function(a){for(var b=!1,c=0;c<d.length;c++)if(angular.equals(f(d[c]),f(a))){b=!0;break}b||d.push(a)}),b=d}return b}}]),angular.module("clotho.interface").directive("contenteditable",["$timeout",function(a){return{require:"?ngModel",link:function(b,c,d,e){function f(){var b=c.text();d.noStripBr||"<br>"!=b||(b=""),b!=e.$modelValue&&e.$setViewValue(b),""===b&&a(function(){c.blur(),c.focus()})}e&&(e.$render=function(){var a=e.$modelValue||"";angular.isObject(a)&&(a=JSON.stringify(a)),c.html(a)},c.on("input",function(){b.$apply(f)}))}}}]),angular.module("clotho.interface").directive("functionCodeDrop",["$upload","$window","$timeout",function(a,b,c){return{restrict:"A",templateUrl:"views/_interface/codeDrop.html",scope:{updateOnRead:"=",showDrop:"@",showBrowser:"@"},controller:["$scope","$element","$attrs",function(a){function d(b,d){b.readAsText(d),b.onload=function(b){c(function(){a.updateOnRead=b.target.result})}}a.uploadRightAway=!1,a.selectedFiles=[],a.inputFiles=[],a.onFileSelect=function(c){a.selectedFiles=c;for(var e=0;e<c.length;e++){var f=c[e];if(b.FileReader&&f.type.indexOf("text")>-1){var g=new FileReader;d(g,c[e])}}}}],link:function(){}}}]),angular.module("clotho.interface").directive("clothoFileUpload",["$upload","$window","$timeout",function(a,b,c){return{restrict:"E",templateUrl:"views/_interface/fileUpload.html",scope:{onSelectCallback:"&?",onSelectEachCallback:"&?",onStart:"&?",startFunction:"&?",uploadUrl:"=?",uploadHeaders:"=?",onUploadSuccess:"&?",onUploadError:"&?"},controller:["$scope","$element","$attrs",function(d,e,f){function g(){d.upload=[],d.uploadResult=[],d.selectedFiles=[],d.started=[],d.dataUrls=[],d.inputFiles=[]}function h(b){return angular.isDefined(d.uploadUrl)?void(d.upload[b]=a.upload({url:d.uploadUrl,headers:d.uploadHeaders,data:{myModel:d.selectedFiles[b].name},file:d.selectedFiles[b],fileFormDataName:"fileUpload"+Date.now().valueOf()}).then(function(a){d.uploadResult.push(a.data.result),angular.isDefined(d.onUploadSuccess)&&d.onUploadSuccess({index:b,result:a.data.result,file:d.selectedFiles[b]})},function(a){console.warn("there was an error uploading",b,d.selectedFiles[b]),angular.isDefined(d.onUploadError)&&d.onUploadError({index:b,result:a,file:d.selectedFiles[b]})})):void console.warn("file could not be uploaded, no url provided")}d.uploadRightAway=!1,d.onFileSelect=function(a){function e(a,b,e){a.onload=function(a){c(function(){d.inputFiles[e]=a.target.result})},angular.isDefined(d.onSelectEachCallback)&&d.onSelectEachCallback({index:e,file:b})}angular.isDefined(d.onSelectCallback)&&d.onSelectEachCallback({files:a});var h;angular.isUndefined(d.selectedFiles)||angular.isDefined(f.replaceOnSelect)?(g(),d.selectedFiles=a,h=0):(h=d.selectedFiles.length,d.selectedFiles=d.selectedFiles.concat(a));for(var i=0;i<a.length;i++,h++){var j=d.selectedFiles[h];if(b.FileReader){var k=new FileReader;j.type.indexOf("text")>-1?(k.readAsText(j),e(k,j,h)):j.type.indexOf("image")>-1?(d.dataUrls[h]=!0,k.readAsDataURL(j),e(k,j,h)):(k.readAsBinaryString(j),e(k,j,h))}d.started[h]=!1,d.uploadRightAway&&d.start(h)}},d.start=function(a){d.started[a]||(d.started[a]=!0,angular.isDefined(d.startFunction)?d.startFunction({index:a,file:d.selectedFiles[a],content:d.inputFiles[a]}):h(a),angular.isDefined(d.onStart)&&d.onStart({index:a,file:d.selectedFiles[a],content:d.inputFiles[a]}))},d.processAll=function(){angular.forEach(d.selectedFiles,function(a,b){d.start(b)})}}],link:function(a,b,c){a.autoUploadVisible=angular.isDefined(c.autoUploadVisible),a.previewFile=function(b){a.showPreview=!0,a.modalFile=a.inputFiles[b],a.modalFileIsImage=!!a.dataUrls[b]}}}}]),angular.module("clotho.interface").directive("backgroundImage",function(){return{restrict:"A",link:function(a,b,c){b.css({"background-position":"50% 50%","background-repeat":"none","background-size":a.$eval(c.backgroundContain)?"contain":"cover"}),c.$observe("ngSrc",function(a){b.css({"background-image":"url("+a+")"})}),a.$watch(function(){return c.backgroundImage},function(a){b.css({"background-image":"url("+a+")"})})}}}),angular.module("clotho.interface").service("$focus",["$document","$timeout","$q","$parse","Clotho","CommandBar",function(a,b,c,d,e,f){var g=f.getCommandBarInput,h=function(b){return Math.max(0,Math.max.apply(null,angular.map(a[0].querySelectorAll(b||"*"),function(a){return parseFloat(angular.element(a).css("z-index"))||null})))},i=function(a,b){return c.when(b.css({"z-index":a,position:"relative"}))},j=function(a){var b=a.css("z-index"),d=h()+1;return i(d,a),c.when(b)},k=function(a,e,f){function g(){h=b(function(){n++,j=e.substring(0,n),i?i.$apply(d(f).assign(i,j)):a[l](j+"|"),g(),n==m&&(i||a[l](a[l]().slice(0,-1)),o.resolve(),b.cancel(h))},Math.round(50*Math.random())+30)}var h,i,j,k={input:!0,textarea:!0,select:!0},l=k[angular.lowercase(a[0].nodeName)]?"val":"text",m=e.length,n=0,o=c.defer();return f&&(i=a.scope()),g(),o.promise},l=function(a){var b=g();return a=angular.isArray(a)?a:[a],b.focus(),k(b,a,f.commandBarInputModel).then(function(){b.remove()})},m=angular.element("<div>").addClass("modal-backdrop fade");m.bind("click",function(a){a.preventDefault(),o()});var n=function(c){return a.find("body").append(m.css("z-index",c||h()+1)),b(function(){m.addClass("in")})},o=function(){return c.when(m.removeClass("in")).then(function(){return b(function(){m.remove()},150)})},p=function(a){var b=a.css("z-index");return n(),i(h()+1,a),c.when(function(){i(b,a),o()
})},q=function(a){var b=c.defer();return a.animate({backgroundColor:"ffff99"},300).animate({backgroundColor:"transparent"},300,function(){b.resolve()}),b.promise};return{maxZ:h,setZ:i,bringToFront:j,typeOut:k,typeOutSearch:l,addBackdrop:n,removeBackdrop:o,highlightElement:p,pulseElement:q}}]),angular.module("clotho.interface").config(["hotkeysProvider",function(a){a.template='<div class="modal fade hotkeyModal" ng-class="{in: helpVisible}"><div class="modal-dialog"><div class="modal-content"><div class="modal-header"><button type="button" class="close" ng-click="toggleCheatSheet()">&times;</button><h4 class="modal-title">{{ title }}</h4></div><div class="modal-body"><table><tbody><tr ng-repeat="hotkey in hotkeys | filter:{ description: \'!$$undefined$$\' }"><td class="hotkey-keys"><span ng-repeat="key in hotkey.format() track by $index" class="hotkey-key">{{ key }}</span></td><td class="hotkey-text">{{ hotkey.description }}</td></tr></tbody></table></div></div><!-- /.modal-content --></div><!-- /.modal-dialog --></div><!-- /.modal -->'}]).run(["hotkeys","$location","$route","Clotho","CommandBar",function(a,b,c,d,e){a.add("f","Focus Command Bar",function(a){a.preventDefault(),e.focusInput()}),a.add("t","Toggle Clotho Terminal",function(){d.trigger("toggleTerminalActive")}),a.add("g h","Go to Homepage",function(){b.path("/")}),c.routes["/browser"]&&a.add("g b","Go to Browser",function(){b.path("/browser")}),c.routes["/editor"]&&a.add("g e","Go to Editor",function(){b.path("/editor")}),c.routes["/trails"]&&a.add("g t","Go to Trails",function(){b.path("/trails")})}]),angular.module("ui.sortable",[]).value("uiSortableConfig",{}).directive("uiSortable",["uiSortableConfig","$log",function(a,b){return{require:"?ngModel",link:function(c,d,e,f){function g(a,b){return b&&"function"==typeof b?function(c,d){a(c,d),b(c,d)}:a}var h={},i={receive:null,remove:null,start:null,stop:null,update:null},j=function(a,b){(b.item.sortable.resort||b.item.sortable.relocate)&&c.$apply()};angular.extend(h,a),f?(f.$render=function(){d.sortable("refresh")},i.start=function(a,b){b.item.sortable={index:b.item.index()}},i.update=function(a,b){b.item.sortable.resort=f},i.receive=function(a,b){b.item.sortable.relocate=!0,f.$modelValue.splice(b.item.index(),0,b.item.sortable.moved)},i.remove=function(a,b){b.item.sortable.moved=1===f.$modelValue.length?f.$modelValue.splice(0,1)[0]:f.$modelValue.splice(b.item.sortable.index,1)[0]},i.stop=function(a,b){if(b.item.sortable.resort&&!b.item.sortable.relocate){var c,d;d=b.item.sortable.index,c=b.item.index(),b.item.sortable.resort.$modelValue.splice(c,0,b.item.sortable.resort.$modelValue.splice(d,1)[0])}},c.$watch(e.uiSortable,function(a){angular.forEach(a,function(a,b){i[b]&&(a=g(i[b],a),"stop"===b&&(a=g(a,j))),d.sortable("option",b,a)})},!0),angular.forEach(i,function(a,b){h[b]=g(a,h[b])}),h.stop=g(h.stop,j)):b.info("ui.sortable: ngModel not provided!",d),d.sortable(h)}}}]),angular.module("clotho.interface").directive("uiMarkdown",function(){var a,b=$clotho.extensions.mixin("bower_components/showdown/src/showdown.js").then(function(){a=new Showdown.converter});return{restrict:"EA",link:function(c,d,e){b.then(function(){if(e.uiMarkdown)c.$watch(e.uiMarkdown,function(b){var c=b?a.makeHtml(b):"";d.html(c)});else{var b=a.makeHtml(d.text());d.html(b)}})}}}),angular.module("clotho.interface").directive("wiki",function(){function a(a){function b(a){return a.replace(/(?:(?:(?:^|\n)[\*#].*)+)/g,function(a){var c=a.match(/(^|\n)#/)?"OL":"UL";return a=a.replace(/(^|\n)[\*#][ ]{0,1}/g,"$1"),a=b(a),"<"+c+"><li>"+a.replace(/^\n/,"").split(/\n/).join("</li><li>")+"</li></"+c+">"})}return console.log(a),a?b(a.replace(/(?:^|\n+)([^# =\*<].+)(?:\n+|$)/gm,function(a,b){return b.match(/^\^+$/)?b:"\n<p>"+b+"</p>\n"}).replace(/(?:^|\n)[ ]{2}(.*)+/g,function(a,b){return b.match(/^\s+$/)?a:"<blockquote>"+b+"</pre>"}).replace(/((?:^|\n)[ ]+.*)+/g,function(a){return a.match(/^\s+$/)?a:"<pre>"+a.replace(/(^|\n)[ ]+/g,"$1")+"</pre>"}).replace(/(?:^|\n)([=]+)(.*)\1/g,function(a,b,c){return"<h"+b.length+">"+c+"</h"+b.length+">"}).replace(/'''(.*?)'''/g,function(a,b){return"<strong>"+b+"</strong>"}).replace(/''(.*?)''/g,function(a,b){return"<em>"+b+"</em>"}).replace(/[^\[](http[^\[\s]*)/g,function(a,b){return'<a href="'+b+'">'+b+"</a>"}).replace(/[\[](http.*)[!\]]/g,function(a,b){var c=b.replace(/[\[\]]/g,"").split(/ /),d=c.shift();return'<a href="'+d+'">'+(c.length?c.join(" "):d)+"</a>"}).replace(/\[\[(.*?)\]\]/g,function(a,b){var c=b.split(/\|/),d=c.shift();return d.match(/^Image:(.*)/)?a:'<a href="'+d+'">'+(c.length?c.join("|"):d)+"</a>"})):a}return{restrict:"EA",link:function(b,c,d){d.wiki?b.$watch(d.wiki,function(b){console.log(b),c.html(a(b))}):c.html(a(c.text()))}}}),angular.module("clotho.interface").controller("SharingCtrl",["$scope","$location","$window","$q","UrlShorten",function(a,b,c,d,e){function f(a){return"images/socialmedia/"+a+".png"}a.social=[{name:"facebook",prefix:"http://www.facebook.com/sharer.php?u=",icon:f("facebook")},{name:"google",prefix:"https://plus.google.com/share?url=",icon:f("google")},{name:"twitter",prefix:"http://twitter.com/share?url=",icon:f("twitter")},{name:"linkedin",prefix:"http://www.linkedin.com/shareArticle?mini=true&url=",icon:f("linkedin")},{name:"digg",prefix:"http://www.digg.com/submit?url=",icon:f("digg")},{name:"reddit",prefix:"http://reddit.com/submit?url=",icon:f("reddit")},{name:"email",prefix:"mailto:?Body=",icon:f("email")}],a.share=function(a,f,g){var h=angular.isDefined(g)?g:b.absUrl(),i=f?e(h):d.when(h);i.then(function(b){var d=a.prefix+b;c.open(d,"email"==a.name?"_self":"_blank")})}}]),angular.module("clotho.urlShorten",[]).factory("UrlShorten",["$q","$http",function(a,b){return function(c){return angular.isUndefined(c)?a.when(null):b.post("https://www.googleapis.com/urlshortener/v1/url",{longUrl:c}).then(function(a){return a.data.id})}}]),angular.module("clotho.interface").directive("darkBackground",["$document",function(a){var b=a.find("body");return{restrict:"A",link:function(a){b.toggleClass("gray",!0),a.$on("$destroy",function(){b.toggleClass("gray",!1)})}}}]);
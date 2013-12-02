angular.module("clotho.interface",["ui.bootstrap","ui.codemirror","ui.keypress"]),angular.module("ui.keypress",[]).factory("keypressHelper",["$parse","$document",function(a,b){var c={8:"backspace",9:"tab",13:"enter",27:"esc",32:"space",33:"pageup",34:"pagedown",35:"end",36:"home",37:"left",38:"up",39:"right",40:"down",45:"insert",46:"delete"},d=function(a){return a.charAt(0).toUpperCase()+a.slice(1)};return function(e,f,g,h){var i,j=[];i=f.$eval(h["ui"+d(e)]),g==b&&(i=h),angular.forEach(i,function(b,c){var d,e;e=a(b),angular.forEach(c.split(" "),function(a){d={expression:e,keys:{}},angular.forEach(a.split("-"),function(a){d.keys[a]=!0}),j.push(d)})});var k=function(a){var b=!(!a.metaKey||a.ctrlKey),d=!!a.altKey,g=!!a.ctrlKey,h=!!a.shiftKey,i=a.keyCode;"keypress"===e&&!h&&i>=97&&122>=i&&(i-=32),angular.forEach(j,function(e){var j=e.keys[c[i]]||e.keys[i.toString()],k=!!e.keys.meta,l=!!e.keys.alt,m=!!e.keys.ctrl,n=!!e.keys.shift;j&&k===b&&l===d&&m===g&&n===h&&f.$apply(function(){e.expression(f,{$event:a})})})};g.bind(e,k),f.$on("$destroy",function(){g.unbind(e,k)})}}]),angular.module("ui.keypress").directive("uiKeydown",["keypressHelper",function(a){return{link:function(b,c,d){a("keydown",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeypress",["keypressHelper",function(a){return{link:function(b,c,d){a("keypress",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeyup",["keypressHelper",function(a){return{link:function(b,c,d){a("keyup",b,c,d)}}}]),angular.module("ui.keypress").service("$keypress",["keypressHelper","$document",function(a,b){return{on:function(c,d,e){return a(c,e,b,d)},off:function(a){a[0].unbind(a[1],a[2])}}}]),angular.module("ui.codemirror",[]).constant("uiCodemirrorConfig",{}).directive("uiCodemirror",["uiCodemirrorConfig",function(a){"use strict";return{restrict:"EA",require:"?ngModel",priority:1,compile:function(b){if(angular.isUndefined(window.CodeMirror))throw new Error("ui-codemirror need CodeMirror to work... (o rly?)");var c=b.text(),d=new CodeMirror(function(a){angular.forEach(b.prop("attributes"),function(b){"ui-codemirror"===b.name?a.setAttribute("ui-codemirror-opts",b.textContent):a.setAttribute(b.name,b.textContent)}),b.parent().length<=0&&b.wrap("<div>"),b.replaceWith(a)},{value:c});return function(b,c,e,f){function g(a){for(var b in a)a.hasOwnProperty(b)&&d.setOption(b,a[b])}var h,i;h=a.codemirror||{},i=angular.extend({},h,b.$eval(e.uiCodemirror),b.$eval(e.uiCodemirrorOpts)),g(i),angular.isDefined(b.$eval(e.uiCodemirror))&&b.$watch(e.uiCodemirror,g,!0),d.on("change",function(a){var c=a.getValue();f&&c!==f.$viewValue&&f.$setViewValue(c),b.$$phase||b.$apply()}),f&&(f.$formatters.push(function(a){if(angular.isUndefined(a)||null===a)return"";if(angular.isObject(a)||angular.isArray(a))throw new Error("ui-codemirror cannot use an object or an array as a model");return a}),f.$render=function(){var a=f.$viewValue||"";d.setValue(a)}),e.uiRefresh&&b.$watch(e.uiRefresh,function(a,b){a!==b&&d.refresh()}),angular.isFunction(i.onLoad)&&i.onLoad(d)}}}}]),angular.module("clotho.interface").provider("$dialog",function(){var a={backdrop:!0,dialogClass:"modal",backdropClass:"modal-backdrop",transitionClass:"fade",triggerClass:"in",resolve:{},backdropFade:!1,dialogFade:!1,keyboard:!0,backdropClick:!0},b={},c={value:0};this.options=function(a){b=a},this.$get=["$http","$document","$compile","$rootScope","$controller","$templateCache","$q","$injector",function(d,e,f,g,h,i,j,k){function l(a){var b=angular.element("<div>");return b.addClass(a),b}function m(c){var d=this,e=this.options=angular.extend({},a,b,c);this._open=!1,this.backdropEl=l(e.backdropClass),e.backdropFade&&(this.backdropEl.addClass(e.transitionClass),this.backdropEl.removeClass(e.triggerClass)),this.modalEl=l(e.dialogClass),e.dialogFade&&(this.modalEl.addClass(e.transitionClass),this.modalEl.removeClass(e.triggerClass)),this.handledEscapeKey=function(a){27===a.which&&(d.close(),a.preventDefault(),d.$scope.$apply())},this.handleBackDropClick=function(a){d.close(),a.preventDefault(),d.$scope.$apply()},this.handleLocationChange=function(){d.close()}}var n=e.find("body");return m.prototype.isOpen=function(){return this._open},m.prototype.open=function(a,b){var c=this,d=this.options;if(a&&(d.templateUrl=a),b&&(d.controller=b),!d.template&&!d.templateUrl)throw new Error("Dialog.open expected template or templateUrl, neither found. Use options or open method to specify them.");return this._loadResolves().then(function(a){var b=a.$scope=c.$scope=a.$scope?a.$scope:g.$new();if(c.modalEl.html(a.$template),c.options.controller){var d=h(c.options.controller,a);c.modalEl.children().data("ngControllerController",d)}f(c.modalEl)(b),c._addElementsToDom(),setTimeout(function(){c.options.dialogFade&&c.modalEl.addClass(c.options.triggerClass),c.options.backdropFade&&c.backdropEl.addClass(c.options.triggerClass)}),c._bindEvents()}),this.deferred=j.defer(),this.deferred.promise},m.prototype.close=function(a){function b(a){a.removeClass(d.options.triggerClass)}function c(){d._open&&d._onCloseComplete(a)}var d=this,e=this._getFadingElements();if(e.length>0){for(var f=e.length-1;f>=0;f--)b(e[f]);return c(),void 0}this._onCloseComplete(a)},m.prototype._getFadingElements=function(){var a=[];return this.options.dialogFade&&a.push(this.modalEl),this.options.backdropFade&&a.push(this.backdropEl),a},m.prototype._bindEvents=function(){this.options.keyboard&&n.bind("keydown",this.handledEscapeKey),this.options.backdrop&&this.options.backdropClick&&this.backdropEl.bind("click",this.handleBackDropClick),this.$scope.$on("$locationChangeSuccess",this.handleLocationChange)},m.prototype._unbindEvents=function(){this.options.keyboard&&n.unbind("keydown",this.handledEscapeKey),this.options.backdrop&&this.options.backdropClick&&this.backdropEl.unbind("click",this.handleBackDropClick)},m.prototype._onCloseComplete=function(a){this._removeElementsFromDom(),this._unbindEvents(),this.$scope.$destroy(),this.deferred.resolve(a)},m.prototype._addElementsToDom=function(){n.append(this.modalEl),this.options.backdrop&&(0===c.value&&n.append(this.backdropEl),c.value++),this._open=!0},m.prototype._removeElementsFromDom=function(){this.modalEl.remove(),this.options.backdrop&&(c.value--,0===c.value&&this.backdropEl.remove()),this._open=!1},m.prototype._loadResolves=function(){var a,b=[],c=[],e=this;if(this.options.dependencies){var f=Application.mixin(this.options.dependencies);c.push("dependencies"),b.push(f)}return this.options.template?a=j.when(this.options.template):this.options.templateUrl&&(a=d.get(this.options.templateUrl,{cache:i}).then(function(a){return a.data})),angular.forEach(this.options.resolve||[],function(a,d){c.push(d),b.push(angular.isString(a)?k.get(a):k.invoke(a))}),c.push("$template"),b.push(a),j.all(b).then(function(a){var b={};return angular.forEach(a,function(a,d){b[c[d]]=a}),b.dialog=e,b})},{dialog:function(a){return new m(a)},messageBox:function(a,b,c){return new m({templateUrl:"views/_interface/ui-custom/dialogMessagebox.html",controller:"MessageBoxController",resolve:{model:function(){return{title:a,message:b,buttons:c}}}})},login:function(){return new m({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,templateUrl:"views/_interface/ui-custom/dialogLogin.html",controller:"DialogLoginController"})},serverAlert:function(a){return new m({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,templateUrl:"views/_interface/ui-custom/dialogMessagebox.html",controller:"ServerAlertController",resolve:{model:function(){return{title:"Server Message",message:a,buttons:[{result:"ok",label:"OK",cssClass:"btn-primary"}]}}}})},share:function(a){return new m({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,templateUrl:"views/_interface/ui-custom/dialogShare.html",controller:"DialogShareController",resolve:{model:function(){return{url:a}}}})},video:function(a,b){return angular.extend(b,{width:"560"}),new m({backdrop:!0,backdropFade:!0,keyboard:!0,backdropClick:!0,template:'<div youtube="{{ videoId }}" params="videoParams""></div>',controller:"VideoDialogController",resolve:{model:function(){return{videoId:a,videoParams:b}}}})}}}]}),angular.module("clotho.interface").controller("MessageBoxController",["$scope","dialog","model",function(a,b,c){a.title=c.title,a.message=c.message,a.buttons=c.buttons,a.close=function(a){b.close(a)}}]),angular.module("clotho.interface").controller("DialogLoginController",["$scope","dialog","Clotho",function(a,b,c){a.close=function(a){b.close(a)},a.notification={},a.cred={username:"",password:""},a.login=function(){c.login(a.cred.username,a.cred.password).then(function(c){console.log("run login"),c?(a.notification={"class":"alert-success",message:"Log in Success"},b.close(a.cred.username)):(a.notification={"class":"alert-error",message:"Log in Error"},a.cred={username:"",password:""})})}}]),angular.module("clotho.interface").controller("ServerAlertController",["$scope","dialog","model","Clotho",function(a,b,c,d){a.title=c.title,a.message=c.message,a.buttons=c.buttons,a.close=function(a){b.close(a)},d.listen("serverAlert",function(){a.close("Another alert appeared"),d.say(a.message)},a)}]),angular.module("clotho.interface").controller("DialogShareController",["$scope","dialog","model","$location",function(a,b,c,d){a.close=function(a){b.close(a)},a.customUrl=c.url&&""!=c.url?c.url:!1,a.social=[{name:"facebook",prefix:"http://www.facebook.com/sharer.php?u="},{name:"google",prefix:"https://plus.google.com/share?url="},{name:"twitter",prefix:"http://twitter.com/share?url="},{name:"linkedin",prefix:"http://www.linkedin.com/shareArticle?mini=true&url="},{name:"digg",prefix:"http://www.digg.com/submit?url="},{name:"reddit",prefix:"http://reddit.com/submit?url="},{name:"email",prefix:"mailto:?Body="}],a.share=function(b){var c=a.customUrl?a.customUrl:b.prefix+d.absUrl();a.close(),window.open(c,"email"==b.name?"_self":"_blank")}}]),angular.module("clotho.interface").controller("VideoDialogController",["$scope","dialog","model",function(a,b,c){a.videoId=c.videoId,a.videoParams=c.videoParams,a.close=function(a){b.close(a)}}]),angular.module("clotho.interface").directive("sharable",["$compile","$http","$templateCache",function(a,b,c){var d=function(a){a.base64icon=base64icon};return{restrict:"E",replace:!0,scope:{content:"="},compile:function(){return{pre:function(d,e){var f;f="undefined"!=typeof d.template?angular.lowercase(d.template):"Instance"!=d.content.type?angular.lowercase(d.content.type):angular.lowercase(/.*\.(.*)$/gi.exec(d.content.schema.name)[1]);var g="views/_interface/sharables/"+f+".html";b.get(g,{cache:c}).error(function(f,h){"404"==h&&b.get("views/_interface/sharables/default.html",{cache:c}).success(function(b){c.put(g,b),e.html(b),a(e.contents())(d)})}).success(function(b){e.html(b),a(e.contents())(d)})},post:d}}}}]),angular.module("clotho.interface").directive("hintButton",function(){return{scope:{hint:"@hintButton"},replace:!0,template:'<button class="btn" popover="{{ hint }}" popover-trigger="mouseenter" popover-placement="left"><i class="icon-info-sign"></i> Hint</button>',link:function(){}}}),angular.module("clotho.interface").directive("closeable",["$compile",function(a){return{restrict:"A",priority:1100,link:function(b,c){b.removeElement=function(){c.remove()},c.prepend(a('<a class="close" style="position: absolute; top: 12px; right: 15px;" ng-click="removeElement()">&times;</a>')(b))}}}]),angular.module("clotho.interface").directive("uiEvent",["$parse",function(a){return function(b,c,d){var e=b.$eval(d.uiEvent);angular.forEach(e,function(d,e){var f=a(d);c.bind(e,function(a){var c=Array.prototype.slice.call(arguments);c=c.splice(1),f(b,{$event:a,$params:c}),b.$$phase||b.$apply()})})}}]),angular.module("clotho.interface").directive("restrictInput",function(){return{restrict:"A",require:"ngModel",link:function(a,b,c,d){var e=function(b){var e=a.$eval(c.restrictInput),f=b.replace(e,"");return f!=b&&(d.$setViewValue(f),d.$render()),f};d.$parsers.unshift(e)}}}),angular.module("clotho.interface").directive("clickOutside",["$document","$parse",function(a,b){return function(c,d,e){var f,g=b(e.clickOutside);e.$observe("clickOutsideActive",function(a){f=!!a});var h=function(a){f&&(a.preventDefault(),a.stopPropagation(),0==d.has(a.target).length&&c.$apply(g(c,{$event:a})))};a.bind("click",h),c.$on("$destroy",function(){a.unbind("click",h)})}}]),angular.module("clotho.interface").filter("breakLines",function(){return function(a,b,c){return(a.match(new RegExp(".{1,"+b+"}","gi"))||[]).join(c||"\n")}}),angular.module("clotho.interface").filter("categorize",["$parse",function(a){return function(b,c){if(c===!1)return b;if((c||angular.isUndefined(c))&&angular.isArray(b)){var d={},e=angular.isString(c)?a(c):function(a){return a},f=function(a){return angular.isObject(a)?e(a):a};angular.forEach(b,function(a){var b=f(a);d[b]||(d[b]=[]),d[b].push(a)}),b=d}return b}}]),angular.module("clotho.interface").filter("highlight",function(){return function(a,b,c,d){return b||angular.isNumber(b)?(c=angular.isDefined(c)?c:"ui-match",a=a.toString(),b=b.toString(),d?a.split(b).join('<span class="'+c+'">'+b+"</span>"):a.replace(new RegExp(b,"gi"),'<span class="'+c+'">$&</span>')):a}}),angular.module("clotho.interface").filter("plainText",function(){return function(a,b){return b?a.replace(/(<([^>]+)>)/gi,""):a}}),angular.module("clotho.interface").filter("shuffle",function(){return function(a,b){if(b===!1)return a;if((b||angular.isUndefined(b))&&angular.isArray(a)){for(var c,d,e=a.slice(0,a.length),f=e.length;f;c=parseInt(Math.random()*f),d=e[--f],e[f]=e[c],e[c]=d);a=e}return a}}),angular.module("clotho.interface").filter("stringEnds",function(){return function(a,b){b=b||15;var c=function(a){return a.substr(0,b)+"..."+a.substring(a.length-b)};if(!angular.isString(a)){if(angular.isArray(a)){var d=angular.isString(a[0])?c(a[0]):"";return"<Array>"+(d?" First Entry:\n"+d:"")}if(angular.isObject(a)){var e=angular.isString(a.key)?c(a.key):"";return"<Object>"+(e?"Object.key:\n"+e:"")}return a}return c(a)}}),angular.module("clotho.interface").filter("unique",["$parse",function(a){return function(b,c){if(c===!1)return b;if((c||angular.isUndefined(c))&&angular.isArray(b)){var d=[],e=angular.isString(c)?a(c):function(a){return a},f=function(a){return angular.isObject(a)?e(a):a};angular.forEach(b,function(a){for(var b=!1,c=0;c<d.length;c++)if(angular.equals(f(d[c]),f(a))){b=!0;break}b||d.push(a)}),b=d}return b}}]),angular.module("clotho.interface").directive("contenteditable",["$timeout",function(a){return{require:"?ngModel",link:function(b,c,d,e){function f(){var b=c.text();d.noStripBr||"<br>"!=b||(b=""),b!=e.$modelValue&&e.$setViewValue(b),""===b&&a(function(){c.blur(),c.focus()})}e&&(e.$render=function(){var a=e.$modelValue||"";angular.isObject(a)&&(a=JSON.stringify(a)),c.html(a)},c.on("input",function(){b.$apply(f)}))}}}]),angular.module("clotho.interface").service("$focus",["$document","$timeout","$q","Clotho",function(a,b,c){var d=$("#clothoCommandBarInput"),e=function(){return Math.max.apply(null,$.map($("body *"),function(a){return"static"!=$(a).css("position")?parseInt($(a).css("z-index"))||1:void 0}))},f=function(a,b){return c.when(b.css({"z-index":a,position:"relative"}))},g=function(a){var b=a.css("z-index"),d=e()+1;return f(d,a),c.when(b)},h=function(a,d,e){function f(a,b,c){for(var d=b.split(".");d.length>1&&(a=a[d.shift()]););a[d.shift()]=c}function g(){h=b(function(){if(l++,a[j](d.substring(0,l)+"|"),g(),l==k){if(a[j](a[j]().slice(0,-1)),e){var c=a.scope();f(c,e,d),c.$apply()}m.resolve(),b.cancel(h)}},Math.round(0*Math.random())+30)}var h,i={input:!0,textarea:!0,select:!0},j=i[angular.lowercase(a[0].nodeName)]?"val":"text",k=d.length,l=0,m=c.defer();return g(),m.promise},i=function(a,e){return c.when(d.focus()).then(function(){return m(d)}).then(function(b){return h(d,a,"display.query").then(function(){return b})}).then(function(a){return b(function(){a()},600).then(function(){return e?(c.when(d.parents("form").submit()),void 0):c.when(d.focus())})})},j=angular.element("<div>").addClass("modal-backdrop fade"),k=function(c){return a.find("body").append(j.css("z-index",c||e()+1)),b(function(){j.addClass("in")})},l=function(){return c.when(j.removeClass("in")).then(function(){return b(function(){j.remove()},150)})},m=function(a){var b=a.css("z-index");return k(),f(e()+1,a),c.when(function(){f(b,a),l()})};return{maxZ:e,setZ:f,bringToFront:g,typeOut:h,typeOutSearch:i,addBackdrop:k,removeBackdrop:l,highlightElement:m}}]),angular.module("ui.sortable",[]).value("uiSortableConfig",{}).directive("uiSortable",["uiSortableConfig","$log",function(a,b){return{require:"?ngModel",link:function(c,d,e,f){function g(a,b){return b&&"function"==typeof b?function(c,d){a(c,d),b(c,d)}:a}var h={},i={receive:null,remove:null,start:null,stop:null,update:null},j=function(a,b){(b.item.sortable.resort||b.item.sortable.relocate)&&c.$apply()};angular.extend(h,a),f?(f.$render=function(){d.sortable("refresh")},i.start=function(a,b){b.item.sortable={index:b.item.index()}},i.update=function(a,b){b.item.sortable.resort=f},i.receive=function(a,b){b.item.sortable.relocate=!0,f.$modelValue.splice(b.item.index(),0,b.item.sortable.moved)},i.remove=function(a,b){b.item.sortable.moved=1===f.$modelValue.length?f.$modelValue.splice(0,1)[0]:f.$modelValue.splice(b.item.sortable.index,1)[0]},i.stop=function(a,b){if(b.item.sortable.resort&&!b.item.sortable.relocate){var c,d;d=b.item.sortable.index,c=b.item.index(),b.item.sortable.resort.$modelValue.splice(c,0,b.item.sortable.resort.$modelValue.splice(d,1)[0])}},c.$watch(e.uiSortable,function(a){angular.forEach(a,function(a,b){i[b]&&(a=g(i[b],a),"stop"===b&&(a=g(a,j))),d.sortable("option",b,a)})},!0),angular.forEach(i,function(a,b){h[b]=g(a,h[b])}),h.stop=g(h.stop,j)):b.info("ui.sortable: ngModel not provided!",d),d.sortable(h)}}}]),$script("bower_components/showdown/src/showdown.js"),angular.module("clotho.interface").directive("uiMarkdown",function(){var a=new Showdown.converter;return{restrict:"EA",link:function(b,c,d){if(d.uiMarkdown)b.$watch(d.uiMarkdown,function(b){var d=b?a.makeHtml(b):"";c.html(d)});else{var e=a.makeHtml(c.text());c.html(e)}}}}),angular.module("clotho.interface").directive("wiki",function(){function a(a){function b(a){return a.replace(/(?:(?:(?:^|\n)[\*#].*)+)/g,function(a){var c=a.match(/(^|\n)#/)?"OL":"UL";return a=a.replace(/(^|\n)[\*#][ ]{0,1}/g,"$1"),a=b(a),"<"+c+"><li>"+a.replace(/^\n/,"").split(/\n/).join("</li><li>")+"</li></"+c+">"})}return console.log(a),a?b(a.replace(/(?:^|\n+)([^# =\*<].+)(?:\n+|$)/gm,function(a,b){return b.match(/^\^+$/)?b:"\n<p>"+b+"</p>\n"}).replace(/(?:^|\n)[ ]{2}(.*)+/g,function(a,b){return b.match(/^\s+$/)?a:"<blockquote>"+b+"</pre>"}).replace(/((?:^|\n)[ ]+.*)+/g,function(a){return a.match(/^\s+$/)?a:"<pre>"+a.replace(/(^|\n)[ ]+/g,"$1")+"</pre>"}).replace(/(?:^|\n)([=]+)(.*)\1/g,function(a,b,c){return"<h"+b.length+">"+c+"</h"+b.length+">"}).replace(/'''(.*?)'''/g,function(a,b){return"<strong>"+b+"</strong>"}).replace(/''(.*?)''/g,function(a,b){return"<em>"+b+"</em>"}).replace(/[^\[](http[^\[\s]*)/g,function(a,b){return'<a href="'+b+'">'+b+"</a>"}).replace(/[\[](http.*)[!\]]/g,function(a,b){var c=b.replace(/[\[\]]/g,"").split(/ /),d=c.shift();return'<a href="'+d+'">'+(c.length?c.join(" "):d)+"</a>"}).replace(/\[\[(.*?)\]\]/g,function(a,b){var c=b.split(/\|/),d=c.shift();return d.match(/^Image:(.*)/)?a:'<a href="'+d+'">'+(c.length?c.join("|"):d)+"</a>"})):a}return{template:"<div></div>",restrict:"EA",scope:{},link:function(b,c,d){d.wiki?b.$watch(d.wiki,function(b,d){d!=b&&c.html(a(b))}):c.html(a(c.text()))}}}),angular.module("clotho.interface").directive("darkBackground",["$document",function(a){return{restrict:"A",link:function(b){var c=a.find("body").css("background");a.find("body").css({background:"#eeeeee"}),b.$on("$destroy",function(){a.find("body").css({background:c})})}}}]);
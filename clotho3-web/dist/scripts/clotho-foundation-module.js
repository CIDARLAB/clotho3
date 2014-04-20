angular.module("clotho.clothoDirectives",["clotho.core","clotho.utils"]),angular.module("clotho.foundation",["clotho.core","clotho.setup","clotho.clothoDirectives","clotho.extensions"]),angular.module("clotho.setup",[]).run(["$rootScope","Clotho",function(a,b){a.Clotho=b}]),angular.module("clotho.utils",["clotho.core"]).service("ClothoUtils",["$q","$http","Clotho",function(a,b,c){function d(a,b){return"widgets/"+a+(b?"/"+b:"")}var e=function(a){return angular.isString(a)&&16==a.length&&/[a-zA-Z0-9]{16}/.test(a)},f=function(a){return b.get(d(a)+"/model.json").then(function(a){return a.data})},g=function(b){var c=[];return _.forEach(b.importedViews,function(a){c.push(f(a).then(function(a){return g(a)}))}),a.all(c).then(function(){var a=[];return _.forEach(b.dependencies,function(c){a.push(d(b.id,c))}),$clotho.extensions.mixin(a)}).then(function(){return b})},h=function(b){function d(a){a.superClass?f.then(function(){return c.get(a.superClass).then(function(a){return e.fields=e.fields.concat(a.fields),d(a)})}):g.resolve()}if(angular.isUndefined(b))return a.when();if(!b.superClass)return a.when(b);var e=angular.copy(b),f=a.when(),g=a.defer();return d(b),g.promise.then(function(){return f}).then(function(){return e})};return{validUUID:e,downloadViewDependencies:g,generateWidgetUrl:d,downloadSchemaDependencies:h}}]),angular.module("clotho.foundation").service("ClothoSchemas",["Clotho","$q",function(a,b){function c(a){return a.schema||null}function d(a){var b=c(a);return b==m||b==n}function e(a){return c(a)==m}function f(a){return c(a)==n}function g(a){return c(a)==p.Function.scaffold.schema}function h(a){return c(a)==p.View.scaffold.schema}function i(a){return d(a)?"Schema":g(a)?"Function":h(a)?"View":"Instance"}function j(a){var b;return b=p[a]?p[a].scaffold:angular.extend({schema:a},p.Instance.scaffold)}function k(a){return angular.map(q,function(b,c){return a[c]})}var l="org.clothocad.core.schema.Schema",m="org.clothocad.core.schema.BuiltInSchema",n="org.clothocad.core.schema.ClothoSchema",o=b.defer(),p={Instance:{readable:"Instance",editor_template_url:"views/_editor/sharable.html",scaffold:{schema:!1,language:"JSONSCHEMA"}},Function:{readable:"Function",editor_template_url:"views/_editor/function.html",scaffold:{schema:"org.clothocad.core.datums.Function",language:"JSONSCHEMA"}},Schema:{readable:"Schema",editor_template_url:"views/_editor/schema.html",scaffold:{schema:n,language:"JSONSCHEMA"}},View:{readable:"View",editor_template_url:"views/_editor/view.html",scaffold:{schema:"org.clothocad.core.datums.View",language:"JSONSCHEMA"}}},q={name:{name:"name",type:"string",description:"Name of the Instance, given by author"},id:{name:"id",type:"string",description:"Unique ID referring to this object"},schema:{name:"schema",type:"string",description:"Pattern describing contents and organization of instance data"},description:{name:"description",type:"string",description:"Description of object, written by author"},author:{name:"author",type:"string",description:"User who created this object"}},r=[{name:"Public",value:"PUBLIC"},{name:"Private",value:"PRIVATE"},{name:"Read Only",value:"READONLY"}],s=[{name:"RegExp",value:"regex"},{name:"Not Null",value:"notnull"}],t={string:"java.lang.String",number:"java.lang.Long","boolean":"java.lang.Boolean",object:"java.util.HashMap",array:"java.util.List"};return a.query({schema:l}).then(function(a){_.remove(a,function(a){return!!p[a.name]}),o.resolve(a)}),{retrievedSchemas:o.promise,sharableTypes:p,accessTypes:r,constraintTypes:s,primitiveToJava:t,isSchema:d,isBuiltIn:e,isClothoSchema:f,determineType:i,determineSchema:c,createScaffold:j,sharableBasicFields:q,pruneToBasicFields:k}}]),angular.module("clotho.clothoDirectives").directive("clothoRun",["Clotho",function(a){var b={input:!0,textarea:!0,select:!0};return{restrict:"A",require:"ngModel",scope:!0,link:function(c,d,e,f){function g(a){k&&f.$setViewValue(a)}function h(a){var b=j?"val":"text";d[b](a)}function i(a){return angular.isArray(a)?a:[a]}var j=!!b[angular.lowercase(d[0].nodeName)];j&&(f.$render=angular.noop);var k=!1;c.$watch(function(){return e.clothoRun},function(a){a&&l(f.$modelValue)}),c.$watch(function(){return f.$modelValue},function(a){l(a)}),c.$watch(function(){return e.clothoRunUpdateModel},function(a){k=!!a});var l=function(b){return b=i(b),a.run(e.clothoRun,b).then(function(a){g(a),h(a)})}}}}]),angular.module("clotho.clothoDirectives").directive("clothoShow",["$q","$http","$timeout","$browser","$rootScope","$compile","Clotho","PubSub","ClothoUtils",function(a,b,c,d,e,f,g,h,i){var j=i.generateWidgetUrl,k=function(a){return b.get(j(a)+"/model.json").then(function(a){return a.data})};return{terminal:!0,restrict:"E",scope:{id:"@",callback:"=?"},controller:["$scope","$element","$attrs",function(){}],link:function(b,e,g){console.log("directive linked"),b.$watch("id",function(l){l&&(console.log(l),e.addClass("clothoWidget"),a.when(k(b.id)).then(function(a){return i.downloadViewDependencies(a)}).then(function(a){if(a.dictionary=angular.extend({},a.dictionary,a.importedViews),a.dictionary.id=a.id,a.bootstrap){var i=a.id+"-additions",k=a.bootstrap.excludeExtensionsModule!==!1?["clotho.extensions"]:[];angular.module(i,k).run(["$rootScope",function(b){angular.extend(b,a.dictionary),b.prefixUrl=function(b,c){return j(c?c:a.id,b)}}]);var l=[];l.push(function(a,b){a.value("$anchorScroll",angular.noop),a.value("$browser",d),a.service("lazyScripts",["$q","$timeout","$document",function(a,b,c){var d=[];this.register=function(a){d.push($clotho.extensions.mixin(a))},b(function(){a.all(d).then(function(){c.triggerHandler("WidgetContentLoaded")})})}]),b.directive("script",["$parse","$rootScope","lazyScripts",function(a,b,c){return{restrict:"E",terminal:!0,compile:function(d,e){if(e.ngSrc){var f=a(e.ngSrc)(b);c.register(f)}}}}])}),l=l.concat(a.bootstrap.modules,i),e.html("<div ng-include=\"prefixUrl('index.html')\"></div>"),a.controller&&g.$set("ng-controller",a.controller),e.data("$injector",null),angular.bootstrap(e,l)}else{var m;if(_.indexOf(a.files,"index.html")>=0)m="<div ",a.controller&&(m+='ng-controller="'+a.controller+'" '),m+="ng-include=\"prefixUrl('index.html')\"></div>";else if(1==a.files.length){var n=/([a-z\-_0-9\/\:\.]*\.(jpg|jpeg|png|gif))/i;n.test(a.files[0])&&(m='<img src="'+j(a.id,a.files[0])+'"alt="view '+a.id+'" />')}angular.extend(b,a.dictionary),e.html(f(m)(b))}c(function(){angular.isFunction(b.callback)&&b.callback(e),h.trigger("clothoShow:"+b.id,[b.id,e,a])})}))})}}}]),angular.module("clotho.clothoDirectives").directive("clothoModal",["$parse","$timeout","hotkeys",function(a,b,c){return{restrict:"E",replace:!0,transclude:!0,templateUrl:"views/_foundation/clothoModal.html",scope:{id:"@?",open:"=?",onClose:"&?",onOpen:"&?",title:"@?"},controller:["$scope","$element","$attrs",function(a){a.$close=function(){a.open&&(a.open=!1,c.del("esc"),b(function(){angular.isFunction(a.onClose)&&a.onClose()}))}}],link:function(a,d,e){e.open||b(function(){a.open=!0}),a.$watch("open",function(b){b&&(a.open=!0,c.add("esc",a.$close),angular.isFunction(a.onOpen)&&a.onOpen())})}}}]),angular.module("clotho.clothoDirectives").directive("sharablePopup",["$animate","$window","$document","$compile","$timeout","Clotho","ClothoSchemas","hotkeys",function(a,b,c,d,e,f,g,h){function i(a){var b=a||"click",c=l[b]||b;return{show:b,hide:c}}var j=c.find("body"),k=function(a){var d=a[0].getBoundingClientRect();return{width:d.width||a.prop("offsetWidth"),height:d.height||a.prop("offsetHeight"),top:d.top+(b.pageYOffset||c[0].body.scrollTop||c[0].documentElement.scrollTop),left:d.left+(b.pageXOffset||c[0].body.scrollLeft||c[0].documentElement.scrollLeft)}},l={mouseenter:"mouseleave",click:"click",focus:"blur"},m="sharablePopup",n='<div sharable-popup-inner sharable="tt_sharable" placement="tt_placement" ></div>';return{restrict:"EA",replace:!0,scope:!0,controller:["$scope","$element","$attrs",function(){}],link:function(b,c,e){function l(){var a,d,e,f;switch(a=k(c),d=v.prop("offsetWidth"),e=v.prop("offsetHeight"),b.tt_placement){case"right":f={top:a.top+a.height/2-e/2,left:a.left+a.width};break;case"top":f={top:a.top-e,left:a.left+a.width/2-d/2};break;case"left":f={top:a.top+a.height/2-e/2,left:a.left-d};break;default:f={top:a.top+a.height,left:a.left+a.width/2-d/2}}f.top+="px",f.left+="px",v.css(f)}function o(){b.tt_isOpen?q():p()}function p(){r()()}function q(){b.$apply(function(){s()})}function r(){return t(),v.css({top:0,left:0,display:"block"}),a.enter(v,j,angular.element(j[0].lastChild),function(){}),l(),b.tt_isOpen=!0,b.$digest(),h.add("esc",q),l}function s(){b.tt_isOpen=!1,u(),h.del("esc")}function t(){v&&u(),v=x(b),b.$digest()}function u(){v&&a.leave(v,function(){v=null})}if(e[m+"Id"]){var v,w,x=d(n),y=!1;c.css({cursor:"pointer"}),b.tt_sharable="loading...",e.$observe(m+"Id",function(a,c){!a||c&&a==c||f.get(a).then(function(a){b.tt_sharable=g.pruneToBasicFields(a)})}),e.$observe(m+"Title",function(a){b.tt_title=a}),e.$observe(m+"Placement",function(a){b.tt_placement=angular.isDefined(a)?a:"bottom"});var z=function(){y&&(c.unbind(w.show,p),c.unbind(w.hide,q))};e.$observe(m+"Trigger",function(a){z(),w=i(a),w.show===w.hide?c.bind(w.show,o):(c.bind(w.show,p),c.bind(w.hide,q)),y=!0}),b.tt_isOpen=angular.isDefined(e[m+"StartOpen"]),b.tt_isOpen,b.$on("$locationChangeSuccess",function(){b.tt_isOpen&&s()}),b.$on("$destroy",function(){z(),u()})}}}}]).directive("sharablePopupInner",function(){return{restrict:"EA",replace:!0,scope:{sharable:"=",placement:"="},templateUrl:"views/_foundation/sharableBasicFieldsPopup.html"}}),!function(a,b,c){function d(a,c){var d=b.createElement("script"),e=k;d.onload=d.onerror=d[p]=function(){d[n]&&!/^c|loade/.test(d[n])||e||(d.onload=d[p]=null,e=1,c())},d.async=1,d.src=a,f.insertBefore(d,f.firstChild)}function e(a,b){q(a,function(a){return!b(a)})}var f=b.getElementsByTagName("head")[0],g={},h={},i={},j={},k=!1,l="push",m="DOMContentLoaded",n="readyState",o="addEventListener",p="onreadystatechange",q=function(a,b){for(var c=0,d=a.length;d>c;++c)if(!b(a[c]))return k;return 1};!b[n]&&b[o]&&(b[o](m,function t(){b.removeEventListener(m,t,k),b[n]="complete"},k),b[n]="loading");var r=function(a,b,f){function k(){if(!--s){g[p]=1,o&&o();for(var a in i)q(a.split("|"),m)&&!e(i[a],m)&&(i[a]=[])}}function m(a){return a.call?a():g[a]}a=a[l]?a:[a];var n=b&&b.call,o=n?b:f,p=n?a.join(""):b,s=a.length;return c(function(){e(a,function(a){j[a]?(p&&(h[p]=1),k()):(j[a]=1,p&&(h[p]=1),d(r.path?r.path+a+".js":a,k))})},0),r};r.get=d,r.ready=function(a,b,c){a=a[l]?a:[a];var d=[];return!e(a,function(a){g[a]||d[l](a)})&&q(a,function(a){return g[a]})?b():!function(a){i[a]=i[a]||[],i[a][l](b),c&&c(d)}(a.join("|")),r};var s=a.$script;r.noConflict=function(){return a.$script=s,this},"undefined"!=typeof module&&module.exports?module.exports=r:a.$script=r}(this,document,setTimeout),angular.module("clotho.extensions",[]).config(["$controllerProvider","$compileProvider","$filterProvider","$provide",function(a,b,c,d){window.$clotho.extensions=$clotho.extensions={},$clotho.extensions.providers={$controllerProvider:a,$compileProvider:b,$filterProvider:c,$provide:d},$clotho.extensions._controller=$clotho.extensions.controller,$clotho.extensions._service=$clotho.extensions.service,$clotho.extensions._factory=$clotho.extensions.factory,$clotho.extensions._value=$clotho.extensions.value,$clotho.extensions._directive=$clotho.extensions.directive,$clotho.extensions.controller=function(b,c){return a.register(b,c),this},$clotho.extensions.service=function(a,b){return d.service(a,b),this},$clotho.extensions.factory=function(a,b){return d.factory(a,b),this},$clotho.extensions.value=function(a,b){return d.value(a,b),this},$clotho.extensions.directive=function(a,c){return b.directive(a,c),this},$clotho.extensions.filter=function(a,b){return c.filter(a,b),this}}]).run(["$rootScope","$q","$timeout","$templateCache","$http","$rootElement","$compile",function(a,b,c,d,e,f){{var g=function(){return angular.module("clotho.extensions")._invokeQueue};g().length}$clotho.extensions.recompile=function(a,b){"undefined"!=typeof a&&(b=b||{},a.hasClass("ng-scope")&&a.scope().$destroy(),f.injector().invoke(function(c,d){var e=d.$new(!0);angular.extend(e,b),c($(a))(e),d.$apply()}))},$clotho.extensions.extend=angular.extend,$clotho.extensions.extendPrimaryRootscope=function(b){$clotho.extensions.extend(a,b)},$clotho.extensions.mixin=function(d){if(angular.isUndefined(d)||""==d)return b.when("no mixin url");var e=b.defer(),f=c(function(){e.reject(null)},5e3);return $script(d,function(){c.cancel(f),a.$safeApply(e.resolve(d))}),e.promise},$clotho.extensions.script=function(a){if(angular.isUndefined(a)||0==a.length)return b.when("no script url");var c=[];return angular.isString(a)&&(a=[a]),angular.forEach(a,function(a){c.push(a+"?_="+Date.now())}),$clotho.extensions.mixin(c)};var h=[];$clotho.extensions.css=function(d){if(angular.isUndefined(d)||""==d)return b.when("no css url");if(_.indexOf(h,d)>-1)return b.when("CSS url already added");var e=b.defer(),f=c(function(){e.reject(null)},5e3);if(document.createStyleSheet)document.createStyleSheet(d),a.$safeApply(e.resolve());else{var g=document.createElement("link");g.type="text/css",g.rel="stylesheet",g.href=d,document.getElementsByTagName("head")[0].appendChild(g),c.cancel(f),a.$safeApply(e.resolve())}return h.push(d),e.promise},$clotho.extensions.cache=function(a){if(angular.isUndefined(a)||""==a)return b.when();var f=b.defer(),g=c(function(){f.reject(null)},5e3);return e.get(a).success(function(b){c.cancel(g),d.put(a,b),f.resolve(b)}).error(function(a){f.reject(a)}),f.promise},$clotho.extensions.bootstrap=angular.bootstrap,$clotho.extensions.determineUrlExtension=function(a){var b=a.split("?")[0];return b.substr(b.lastIndexOf(".")+1)};document.getElementsByTagName("head")[0]}]);
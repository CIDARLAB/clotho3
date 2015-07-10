angular.module("clotho.clothoDirectives",["clotho.core","clotho.utils"]),angular.module("clotho.foundation",["clotho.core","clotho.setup","clotho.clothoDirectives","clotho.extensions"]),angular.module("clotho.setup",[]).run(["$rootScope","Clotho",function(a,b){a.Clotho=b}]),angular.module("clotho.utils",["clotho.core"]).service("ClothoUtils",["$q","$http","Clotho",function(a,b){function c(a,b){return"widgets/"+a+(b?"/"+b:"")}var d=function(a){return angular.isString(a)&&16==a.length&&/[a-zA-Z0-9]{16}/.test(a)},e=function(a){return b.get(c(a)+"/model.json").then(function(a){return a.data})},f=e,g=function(b){var d=[];return angular.forEach(b.importedViews,function(a){d.push(f(a).then(function(a){return g(a)}))}),a.all(d).then(function(){var a=[];return angular.forEach(b.dependencies,function(d){a.push(c(b.id,d))}),$clotho.extensions.mixin(a)}).then(function(){return b})};return{validUUID:d,getViewInfo:f,downloadViewDependencies:g,generateWidgetUrl:c}}]),angular.module("clotho.foundation").service("ClothoSchemas",["Clotho","Debug","$filter","$q",function(a,b,c,d){function e(a){return null===a?"null":angular.isBoolean(a)?"boolean":angular.isNumber(a)?"number":angular.isArray(a)?"array":angular.isObject(a)?"object":"string"}function f(a){return B[a]?B[a].class:"default"}function g(a){return A[a]||""}function h(b){var e=b.superClass;return e?z.promise.then(function(b){var f=c("filter")(b,{id:e})[0];return angular.isEmpty(f)?d.when(f):a.get(e)}):d.when(null)}function i(a){return angular.isObject(a)&&!angular.isEmpty(a)&&angular.isDefined(a.id)&&angular.isDefined(a.schema)}function j(a){return!angular.isEmpty(a)&&angular.isDefined(a.schema)?a.schema:null}function k(a){var b=j(a);return b==x||b==y||b==w}function l(a){return j(a)==x}function m(a){return j(a)==y}function n(a){var b=j(a);return b==B.Function.schema||"org.clothocad.core.datums.Module"==b}function o(a){return j(a)==B.View.schema}function p(b){return a.run("clotho.functions.schema.determineSharableType",[b.id],{mute:!0})}function q(a){return k(a)?"Schema":n(a)?"Function":o(a)?"View":i(a)?"Instance":""}function r(a,b){return angular.isUndefined(b)?!1:"org.clothocad.model.NucSeq"==b?/clotho.demo.sequence.*/.test(a.id)||_.indexOf(["org.clothocad.model.NucSeq","Part","Vector"],j(a))>0:j(a)==b}function s(a){return M[a]||M["default"]}function t(a){var b;return b=B[a]?angular.extend({schema:B[a].schema},B[a].scaffold):angular.extend({schema:a},B.Instance.scaffold)}function u(a){return angular.map(C,function(b,c){return a[c]})}var v=new b("clothoSchemas","#992299"),w="org.clothocad.core.schema.Schema",x="org.clothocad.core.schema.BuiltInSchema",y="org.clothocad.core.schema.ClothoSchema",z=d.defer(),A={},B={Instance:{readable:"Instance",editor_template_url:"views/_editor/sharable.html",schema:!1,scaffold:{language:"JSONSCHEMA"},"class":"info"},Function:{readable:"Function",editor_template_url:"views/_editor/function.html",schema:"org.clothocad.core.datums.Function",scaffold:{},"class":"success"},Schema:{readable:"Schema",editor_template_url:"views/_editor/schema.html",schema:y,scaffold:{language:"JSONSCHEMA"},"class":"danger"},View:{readable:"View",editor_template_url:"views/_editor/view.html",schema:"org.clothocad.core.datums.View",scaffold:{},"class":"warning"}},C={name:{name:"name",type:"string",description:"Name of the Instance, given by author"},id:{name:"id",type:"string",description:"Unique ID referring to this object"},schema:{name:"schema",type:"string",description:"Pattern describing contents and organization of instance data"},description:{name:"description",type:"string",description:"Description of object, written by author"},author:{name:"author",type:"string",description:"User who created this object"}},D=[{name:"Public",value:"PUBLIC"},{name:"Private",value:"PRIVATE"},{name:"Read Only",value:"READONLY"}],E=[{name:"RegExp",value:"regex"},{name:"Not Null",value:"notnull"}],F={"boolean":"java.lang.Boolean",number:"java.lang.Long",object:"java.util.HashMap",array:"java.util.List",string:"java.lang.String"},G=/[^\s"']+|"([^"]*)"|'([^']*)'/,H=(G.toString(),{"boolean":{type:"boolean",input:{type:"checkbox"}},number:{type:"number",input:{type:"number"}},string:{type:"string",input:{type:"text"}},date:{type:"date",input:{type:"date"}},array:{type:"array"},"org.bson.types.ObjectId":{type:"id",input:{type:"text"}},object:{type:"object"},"default":{type:"string",input:{type:"text"}}}),I={"java.lang.Boolean":{type:"boolean",input:{type:"checkbox"}},"java.lang.Number":{type:"number",input:{type:"number"}},"java.lang.Long":{type:"number",input:{type:"number"}},"java.lang.Integer":{type:"number",input:{type:"number","ng-pattern":"/^[-+]?\\d+$/"}},"java.lang.Short":{type:"number",input:{type:"number",min:-32767,max:32767}},"java.lang.String":{type:"string",input:{type:"text"}},"java.awt.Color":{type:"color",input:{type:"color"}},"java.util.Date":{type:"date",input:{type:"date"}},"java.util.Set":{type:"array"},"java.util.List":{type:"array"},"org.bson.types.ObjectId":{type:"id",input:{type:"text"}},"java.util.Map":{type:"object"},"java.util.HashMap":{type:"object"},"default":{type:"string",input:{type:"text"}}};a.query({schema:w},{mute:!0}).then(function(a){angular.forEach(a,function(a){A[a.id]=a.name}),z.resolve(a)});var J=function(b){return a.run("clotho.functions.schema.getParents",[b],{mute:!0})},K=function(a){function b(a){a.superClass?e.then(function(){return h(a).then(function(a){return c=c.concat(a.fields),b(a)},function(){v.warn("couldnt get parent schema "+a.superClass),f.reject(c)})}):f.resolve()}var c=[];if(angular.isUndefined(a)||!a.superClass)return d.when(c);var e=d.when(),f=d.defer();return b(a),f.promise.then(function(){return e}).then(function(){return c})},L=function(b){function c(b){b.superClass?f.then(function(){return a.get(b.superClass).then(function(a){return e.fields=e.fields.concat(a.fields),c(a)},function(){v.warn("couldnt get parent schema "+b.superClass),g.reject(e)})}):g.resolve()}if(angular.isUndefined(b))return d.when();if(!b.superClass)return d.when(b);var e=angular.copy(b),f=d.when(),g=d.defer();return c(b),g.promise.then(function(){return f}).then(function(){return e})},M={Instance:"glyphicon glyphicon-file",Function:"glyphicon glyphicon-play-circle",View:"glyphicon glyphicon-picture",Schema:"glyphicon glyphicon-cog","default":"glyphicon glyphicon-file"};return{retrievedSchemas:z.promise,downloadSchemaDependencies:L,getSuperclassFields:K,getParentSchemaIds:J,sharableTypes:B,isBasicField:function(a){return angular.isDefined(C[angular.lowercase(a)])},isPrimitiveField:function(a){return angular.isDefined(F[angular.lowercase(a)])},accessTypes:D,constraintTypes:E,primitiveToJava:F,formTypeMap:H,javaToJavascript:I,isSharable:i,isFunction:n,isView:o,isSchema:k,isInstance:function(a){return"Instance"==q(a)},isBuiltIn:l,isClothoSchema:m,mapSchemaIdToName:g,isInstanceOfSchema:r,dirtyDetermineType:q,determineSharableType:p,determineSharableIcon:s,determineFieldType:e,determineSchema:j,createScaffold:t,typeToColorClass:f,sharableBasicFields:C,pruneToBasicFields:u}}]),angular.module("clotho.clothoDirectives").directive("clothoRun",["Clotho","$parse",function(a,b){var c={input:!0,textarea:!0,select:!0};return{restrict:"A",require:"ngModel",scope:!0,link:function(d,e,f,g){function h(a){angular.isDefined(l.assign)&&l.assign(d,a)}function i(a){var b=k?"val":"text";e[b](a)}function j(a){return angular.isArray(a)?a:[a]}var k=!!c[angular.lowercase(e[0].nodeName)];k&&(g.$render=angular.noop);var l=b(f.clothoRunUpdate);d.$watch(function(){return f.clothoRun},function(a){a&&m(g.$modelValue)}),d.$watch(function(){return g.$modelValue},function(a){m(a)});var m=function(b){return b=j(b),a.run(f.clothoRun,b).then(function(a){h(a),i(a)})}}}}]),angular.module("clotho.clothoDirectives").directive("clothoShow",["$q","$http","$timeout","$browser","$rootScope","$compile","Clotho","PubSub","ClothoUtils",function(a,b,c,d,e,f,g,h,i){var j=i.generateWidgetUrl,k=/([a-z\-_0-9\/\:\.]*\.(jpg|jpeg|png|gif))/i;return{terminal:!0,restrict:"E",scope:{id:"@",callback:"=?"},controller:["$scope","$element","$attrs",function(){}],link:function(b,e,g){b.$watch("id",function(l){l&&(e.addClass("clothoWidget"),a.when(i.getViewInfo(b.id)).then(function(a){return console.log(a),i.downloadViewDependencies(a)}).then(function(a){if(a.dictionary=angular.extend({},a.dictionary,a.importedViews,{id:a.id}),a.bootstrap){var i=a.id+"-additions",l=a.bootstrap.excludeExtensionsModule!==!1?["clotho.extensions"]:[];angular.module(i,l).run(["$rootScope",function(b){angular.extend(b,a.dictionary),b.prefixUrl=function(b,c){return j(c?c:a.id,b)}}]);var m=[];m.push(function(a,b){a.value("$anchorScroll",angular.noop),a.value("$browser",d),a.service("lazyScripts",["$q","$timeout","$document",function(a,b,c){var d=[];this.register=function(a){d.push($clotho.extensions.mixin(a))},b(function(){a.all(d).then(function(){c.triggerHandler("WidgetContentLoaded")})})}]),b.directive("script",["$parse","$rootScope","lazyScripts",function(a,b,c){return{restrict:"E",terminal:!0,compile:function(d,e){if(e.ngSrc){var f=a(e.ngSrc)(b);c.register(f)}}}}])}),m=m.concat(a.bootstrap.modules,i),e.html("<div ng-include=\"prefixUrl('index.html')\"></div>"),a.controller&&g.$set("ng-controller",a.controller),e.data("$injector",null),angular.bootstrap(e,m)}else{b.prefixUrl=function(b,c){return j(c?c:a.id,b)};var n;_.indexOf(a.files,"index.html")>=0?(n="<div ",a.controller&&(n+='ng-controller="'+a.controller+'" '),n+="ng-include=\"prefixUrl('index.html')\"></div>"):1==a.files.length&&k.test(a.files[0])&&(n='<img src="'+j(a.id,a.files[0])+'"alt="view '+a.id+'" />'),angular.extend(b,a.dictionary),e.html(f(n)(b))}c(function(){angular.isFunction(b.callback)&&b.callback(e),h.trigger("clothoShow:"+b.id,[b.id,e,a])})}))})}}}]),angular.module("clotho.clothoDirectives").directive("clothoModal",["$parse","$timeout","$http","$compile","$sce","$injector",function(a,b,c,d,e,f){var g,h=f.has("hotkeys");return h&&(g=f.get("hotkeys")),{restrict:"E",replace:!0,transclude:!0,templateUrl:"views/_foundation/clothoModal.html",scope:{id:"@?",open:"=?",onClose:"&?",onOpen:"&?",title:"@?",content:"=?",templateUrl:"=?",actions:"=?",modalSize:"@"},controller:["$scope","$element","$attrs",function(a,c){a.$close=function(){a.open&&(a.open=!1,g.del("esc"),b(function(){a.onClose()}))},a.clickCheckBackdrop=function(b){b.target==c[0]&&a.$close()}}],link:function(a,f,i){i.open||b(function(){a.open=!0}),a.$watch("content",function(b){a.contentTrusted=e.trustAsHtml(b)}),a.$watch("templateUrl",function(b,e){!b||b==e&&e&&a.hasTemplate||c.get(b,{cache:!0}).success(function(b){angular.element(f[0].querySelector("[template-insert]")).html(d(b)(a)),a.hasTemplate=!0}).error(function(){a.hasTemplate=!1})}),a.$watch("open",function(b){b&&(a.open=!0,h&&g.bindTo(a).add("esc",a.$close),angular.isFunction(a.onOpen)&&a.onOpen())}),a.$on("$destroy",function(){b(function(){a.onClose()})})}}}]).service("$clothoModal",["$window","$rootScope","$compile","$parse",function(a,b,c,d){function e(){i(h),h&&h.$destroy(),g&&g.remove(),g=null,h=null,i=angular.noop}var f=angular.element(a.document.body),g=null,h=null,i=angular.noop;this.create=function(a,i){e();var j=angular.extend({},a);h=angular.isDefined(i)?i.$new():b.$new();var k=d(j["on-close"])||angular.noop;h.clothoModalClose=function(){k(h),e()},j["on-close"]="clothoModalClose()",g=c(angular.element("<clotho-modal>").attr(j))(h),f.append(g)},this.destroy=e}]),angular.module("clotho.clothoDirectives").provider("$clothoPopup",function(){function a(a){var b=a||"click",d=c[b]||b;return{show:b,hide:d}}function b(a){var b=/[A-Z]/g,c="-";return a.replace(b,function(a,b){return(b?c:"")+a.toLowerCase()})}var c={mouseenter:"mouseleave",click:"click",focus:"blur"},d={};this.$get=["$animate","$injector","$window","$document","$compile","$timeout","$parse","Clotho","ClothoSchemas",function(c,e,f,g,h,i){var j,k=e.has("hotkeys");return k&&(j=e.get("hotkeys")),function(e,l){l=angular.extend({},d,l);var m=g.find("body"),n=function(a){var b=a[0].getBoundingClientRect();return{width:b.width||a.prop("offsetWidth"),height:b.height||a.prop("offsetHeight"),top:b.top+(f.pageYOffset||g[0].documentElement.scrollTop),left:b.left+(f.pageXOffset||g[0].documentElement.scrollLeft)}},o="<div "+b(e)+'-inner sharable-id="sharable_id" sharable-model="passedModel" popup-title="popup_title" popup-placement="{{popup_placement}}" reposition="repositionFunction()"></div>';return{restrict:"EA",scope:{passedModel:"=?"+e+"Model"},compile:function(){var b=h(o);return function(d,f,g){function h(){var a,b,c,e;switch(a=n(f),b=x.prop("offsetWidth"),c=x.prop("offsetHeight"),d.popup_placement){case"topRight":e={top:a.top+a.height/2-55,left:a.left+a.width};break;case"topLeft":e={top:a.top+a.height/2-55,left:a.left-b};break;case"right":e={top:a.top+a.height/2-c/2,left:a.left+a.width};break;case"top":e={top:a.top-c,left:a.left+a.width/2-b/2};break;case"left":e={top:a.top+a.height/2-c/2,left:a.left-b};break;case"bottom":e={top:a.top+a.height,left:a.left+a.width/2-b/2};break;default:e={top:a.top+a.height,left:a.left}}e.top+="px",e.left+="px",x.toggleClass("right","topRight"==d.popup_placement),x.toggleClass("left","topLeft"==d.popup_placement),x.css(e)}function o(){return d.popupForceOpen||d.popupOpen}function p(){d.popupOpen?r():q()}function q(){t()()}function r(){d.$apply(function(){u()})}function s(a){angular.isDefined(a)?a?t()():u():o()?t()():u()}function t(){return angular.isEmpty(d.sharable_id)&&angular.isEmpty(d.passedModel)?(d.popupOpen=!1,angular.noop):(v(),x.css({top:0,left:0,display:"block"}),c.enter(x,m,angular.element(m[0].lastChild),angular.noop),h(),d.popupOpen=!0,k&&j.bindTo(y).add("esc",u),h)}function u(a){(a||!d.popupForceOpen)&&(d.popupOpen=!1,x&&w())}function v(){x&&w(),y=d.$new(),x=b(y,function(){}),y.$digest()}function w(){x&&(x.remove(),x=null),y&&(y.$destroy(),y=null)}if(!angular.isUndefined(g[e+"Id"])||!angular.isUndefined(g[e+"Model"])){var x,y,z,A=!1;d.repositionFunction=function(){i(function(){h()})},g.$observe(e+"Id",function(a,b){!a||b&&a==b||(d.sharable_id=a)}),g.$observe(e+"Title",function(a,b){!a||b&&a==b||(d.popup_title=a)}),g.$observe(e+"Placement",function(a){d.popup_placement=a||l.placement});var B=function(){A&&(f.unbind(z.show,q),f.unbind(z.hide,r))};g.$observe(e+"Trigger",function(b){B(),"none"!=b&&(z=a(b||l.trigger),z.show===z.hide?f.bind(z.show,p):(f.bind(z.show,q),f.bind(z.hide,r)),A=!0)}),d.$watch(function(){return g[e+"Open"]},function(a){d.popupForceOpen=d.$eval(a),d.popupOpen=d.$eval(a),setTimeout(function(){s(d.popupOpen)})}),d.$on("$locationChangeSuccess",function(){d.popupOpen&&u(!0)}),d.$on("$destroy",function(){B(),u(!0)})}}}}}}]}),angular.module("clotho.clothoDirectives").directive("sharablePopup",["$clothoPopup",function(a){return a("sharablePopup")}]).directive("sharablePopupInner",["Clotho","ClothoSchemas","$timeout","$injector",function(a,b,c,d){var e=d.has("clothoEditorDirective");return{restrict:"EA",replace:!0,scope:{sharableId:"=?",sharableModel:"=?",placement:"@popupPlacement",reposition:"&"},templateUrl:"views/_foundation/sharableBasicFieldsPopup.html",link:function(d){function f(c){d.sharable=c,b.determineSharableType(c).then(function(a){d.type=a},function(){d.type=b.dirtyDetermineType(c)}).then(function(){d.iconClass=b.determineSharableIcon(d.type),d.labelClass="label-"+b.typeToColorClass(d.type)}),b.isSchema(c)&&(d.isSchema=!0,a.get(c.id,{mute:!0}).then(function(a){b.getSuperclassFields(a).then(function(b){d.schema=a,d.inheritedFields=b})}))}d.editorPresent=e,d.$watch("sharableModel",function(b){b&&(f(b),b.id&&a.get(b.id,{mute:!0}).then(function(a){d.fullSharable=a}))}),d.$watch("sharableId",function(c){c&&a.get(c,{mute:!0}).then(function(a){angular.isEmpty(a)||(d.fullSharable=a,f(b.pruneToBasicFields(a)),d.reposition())})}),d.showView=function(){d.activeView=d.activeView?"":d.type,c(function(){d.reposition()})},d.edit=a.edit,d.$on("$destroy",function(){})}}}]),!function(a,b,c){function d(a,c){var d=b.createElement("script"),e=k;d.onload=d.onerror=d[p]=function(){d[n]&&!/^c|loade/.test(d[n])||e||(d.onload=d[p]=null,e=1,c())},d.async=1,d.src=a,f.insertBefore(d,f.firstChild)}function e(a,b){q(a,function(a){return!b(a)})}var f=b.getElementsByTagName("head")[0],g={},h={},i={},j={},k=!1,l="push",m="DOMContentLoaded",n="readyState",o="addEventListener",p="onreadystatechange",q=function(a,b){for(var c=0,d=a.length;d>c;++c)if(!b(a[c]))return k;return 1};!b[n]&&b[o]&&(b[o](m,function t(){b.removeEventListener(m,t,k),b[n]="complete"},k),b[n]="loading");var r=function(a,b,f){function k(){if(!--s){g[p]=1,o&&o();for(var a in i)q(a.split("|"),m)&&!e(i[a],m)&&(i[a]=[])}}function m(a){return a.call?a():g[a]}a=a[l]?a:[a];var n=b&&b.call,o=n?b:f,p=n?a.join(""):b,s=a.length;return c(function(){e(a,function(a){j[a]?(p&&(h[p]=1),k()):(j[a]=1,p&&(h[p]=1),d(r.path?r.path+a+".js":a,k))})},0),r};r.get=d,r.ready=function(a,b,c){a=a[l]?a:[a];var d=[];return!e(a,function(a){g[a]||d[l](a)})&&q(a,function(a){return g[a]})?b():!function(a){i[a]=i[a]||[],i[a][l](b),c&&c(d)}(a.join("|")),r};var s=a.$script;r.noConflict=function(){return a.$script=s,this},"undefined"!=typeof module&&module.exports?module.exports=r:a.$script=r}(this,document,setTimeout),$clotho.extensions=angular.module("clotho.extensions",[]).config(["$controllerProvider","$compileProvider","$filterProvider","$provide","$injector","$animateProvider",function(a,b,c,d,e,f){$clotho.extensions.providers={$controllerProvider:a,$compileProvider:b,$filterProvider:c,$provide:d,$injector:e,$animateProvider:f},$clotho.extensions.controller=function(){return a.register.apply(this,arguments),this},$clotho.extensions.service=function(){return d.service.apply(this,arguments),this},$clotho.extensions.factory=function(){return d.factory.apply(this,arguments),this},$clotho.extensions.value=function(){return d.value.apply(this,arguments),this},$clotho.extensions.directive=function(){return b.directive.apply(this,arguments),this},$clotho.extensions.filter=function(){return c.register.apply(this,arguments),this},$clotho.extensions.animation=function(){return f.register.apply(this,arguments),this}}]).run(["ClothoExtensions",function(a){angular.isFunction(window.$clotho.extensions.downloadDependencies)||angular.extend(window.$clotho.extensions,a)}]).service("ClothoExtensions",["$rootScope","$q","$timeout","$templateCache","$http","$cacheFactory","$rootElement",function(a,b,c,d,e,f,g){function h(){return angular.module("clotho.extensions")._invokeQueue}function i(){s=h().length}function j(a){var b=(new Date).getTime();return a.indexOf("?")>=0?"&"===a.substring(0,a.length-1)?a+"_dc="+b:a+"&_dc="+b:a+"?_dc="+b}function k(a){var b=a.split("?")[0];return b.substr(b.lastIndexOf(".")+1)}function l(a,b){$script(a,function(){b(a)})}function m(a,b){l(j(a),b)}function n(a,b){if(document.createStyleSheet)document.createStyleSheet(a);else{var c=document.createElement("link");c.type="text/css",c.rel="stylesheet",c.href=a,document.getElementsByTagName("head")[0].appendChild(c)}b(a)}function o(a,b,c){e.get(a,b).success(function(e){d.put(b.forceName||a,e),c(a)}).error(function(){c(null)})}function p(d,e){function f(){c.cancel(h),e.noCache!==!1&&u.put(e.cacheName||d),a.$evalAsync(g.resolve(d))}if(angular.isUndefined(d)||0==d.length)return b.when(null);if(e=angular.extend({mixin:!0},e),e.mixin!==!1&&u.get(d))return b.when(d);var g=b.defer(),h=c(function(){g.reject(null)},t);switch(e.forceType||k(d)){case"js":e.mixin?l(d,f):m(d,f);break;case"css":n(d,f);break;case"html":o(d,e,f);break;default:g.reject(new Error("Dont know how to handle "+d))}return g.promise}function q(a,c){a=angular.isArray(a)?a:[a];var d=[];return angular.forEach(a,function(a){d.push(p(a,c))}),b.all(d)}var r={},s=0,t=5e3,u=f("fileUrls");return i(),r.downloadDependencies=function(a){return angular.isObject(a)&&Object.keys(a).length>0?b.all([r.css(a.css),r.mixin(a.mixin),r.script(a.script)]).then(function(){return function(){c(function(){r.script(a.onload)})}}):q(a)},r.mixin=q,r.script=function(a,b){return q(a,angular.extend({},b,{mixin:!1}))},r.css=q,r.cache=function(a,b){return b=angular.isString(b)?{templateName:b}:angular.extend({},b),q(a,b)},r.determineUrlExtension=k,r.recompile=function(a,b){angular.isUndefined(a)||(b=b||{},a=angular.element(a),a.scope()&&a.scope().$destroy(),g.injector().invoke(["$compile","$rootScope",function(c,d){var e=d.$new(!0);angular.extend(e,b),c(a)(e),d.$apply()}]))},r.extendRootscope=function(b){a.$evalAsync(angular.extend(a,b))},r}]),angular.module("clotho.foundation").directive("clothoShowAuth",["PubSub","ClothoAuth",function(a,b){function c(a,b){var c=a.$eval(b);return!angular.isString(c)&&angular.isArray(c)&&(c=b),angular.isString(c)&&(c=c.split(",")),c}function d(a,b){for(var c=0;c<b.length;c++)if(b[c]===a)return!0;return!1}function e(a){if(!a.length)throw new Error("clotho-show-auth directive must be (you may use a comma-separated list): "+f.join(", "));return angular.forEach(a,function(a){if(!d(a,f))throw new Error('Invalid state(s) "'+a+'" for clotho-show-auth directive, must be one of: '+f.join(", "))}),!0}var f=["login","logout","error"];return{restrict:"A",link:function(a,f,g){function h(a){var b=d(a,i);setTimeout(function(){f.toggleClass("ng-cloak",!b)},0)}var i=c(a,g.clothoShowAuth);e(i),b.addStateListener(h)}}}]),angular.module("clotho.foundation").service("Facebook",["$q","$window","$filter",function(a,b,c){function d(){g={},h=a.defer(),i=h.promise}function e(){l=!0;var b=a.defer();return o.then(function(){FB.api("/me",function(a){f(a),b.resolve(a)})}),b.promise}function f(a){l=!0,m=!1,angular.extend(g,a,{provider:"facebook",icon:"//graph.facebook.com/"+a.id+"/picture?width=500&height=500"}),h.resolve(g)}var g,h,i,j=this,k=0x52d209de19bff,l=!1,m=!1,n=a.defer(),o=n.promise;d(),b.fbAsyncInit=function(){FB.init({appId:k,xfbml:!0,version:"v2.0",cookie:!0,status:!0}),FB.getLoginStatus(function(a){"connected"===a.status?e():"not_authorized"===a.status,n.resolve(FB)}),FB.Event.subscribe("auth.login",function(){}),FB.Event.subscribe("auth.logout",function(){})},j.login=function(b){var c=angular.extend({scope:"public_profile email"},b),d=a.defer();return o.then(function(a){l?i.then(function(a){d.resolve(a)}):a.login(function(a){console.log("FB login response",a),a.authResponse?e().then(function(){i.then(function(a){d.resolve(a)})}):d.reject(a)},c)}),d.promise},j.logout=function(a){return o.then(function(b){b.getLoginStatus(function(c){"connected"===c.status&&b.logout(a)}),d()})},j.getUser=function(){var b=a.defer();return angular.isEmpty(g)?o.then(function(){l?i.then(function(a){b.resolve(a)}):b.reject()}):b.resolve(g),b.promise},j.convertToPersonSharable=function(a){return a=a||g,angular.isEmpty(g)?void 0:{schema:"org.clothocad.model.Person",id:a.email,name:a.name,emailAddress:a.email,dateCreated:c("date")(Date.now().valueOf(),"yyyy-MM-dd"),icon:a.icon,social:{facebook:a.link}}},function(a,b,c){var d,e=a.getElementsByTagName(b)[0];a.getElementById(c)||(d=a.createElement(b),d.id=c,d.src="//connect.facebook.net/en_US/sdk.js",e.parentNode.insertBefore(d,e))}(document,"script","facebook-jssdk")}]);
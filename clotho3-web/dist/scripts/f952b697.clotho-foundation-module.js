angular.module("clotho.clothoDirectives",["clotho.core","clotho.utils"]),angular.module("clotho.foundation",["clotho.core","clotho.setup","clotho.clothoDirectives","clotho.extensions"]),angular.module("clotho.setup",[]).run(["$rootScope","Clotho",function(a,b){a.Clotho=b}]),angular.module("clotho.utils",["clotho.core"]).service("ClothoUtils",["$q","$http","Clotho",function(a,b){function c(a,b){return"widgets/"+a+(b?"/"+b:"")}var d=function(a){return angular.isString(a)&&16==a.length&&/[a-zA-Z0-9]{16}/.test(a)},e=function(a){return b.get(c(a)+"/model.json").then(function(a){return a.data})},f=function(b){var d=[];return _.forEach(b.importedViews,function(a){d.push(e(a).then(function(a){return f(a)}))}),a.all(d).then(function(){var a=[];return _.forEach(b.dependencies,function(d){a.push(c(b.id,d))}),$clotho.extensions.mixin(a)}).then(function(){return b})};return{validUUID:d,downloadViewDependencies:f,generateWidgetUrl:c}}]),angular.module("clotho.foundation").service("ClothoSchemas",["Clotho","Debug","$filter","$q",function(a,b,c,d){function e(a){return null===a?"null":angular.isBoolean(a)?"boolean":angular.isNumber(a)?"number":angular.isArray(a)?"array":angular.isObject(a)?"object":"string"}function f(a){return B[a]?B[a].class:"default"}function g(a){return A[a]||""}function h(b){var e=b.superClass;return e?z.promise.then(function(b){var f=c("filter")(b,{id:e})[0];return angular.isEmpty(f)?d.when(f):a.get(e)}):d.when(null)}function i(a){return angular.isObject(a)&&!angular.isEmpty(a)&&angular.isDefined(a.id)&&angular.isDefined(a.schema)}function j(a){return!angular.isEmpty(a)&&angular.isDefined(a.schema)?a.schema:null}function k(a){var b=j(a);return b==x||b==y}function l(a){return j(a)==x}function m(a){return j(a)==y}function n(a){var b=j(a);return b==B.Function.schema||"org.clothocad.core.datums.Module"==b}function o(a){return j(a)==B.View.schema}function p(b){return a.run("clotho.functions.schema.determineSharableType",[b.id],{mute:!0})}function q(a){return k(a)?"Schema":n(a)?"Function":o(a)?"View":i(a)?"Instance":""}function r(a,b){return angular.isUndefined(b)?!1:"org.clothocad.model.NucSeq"==b?/clotho.demo.sequence.*/.test(a.id)||_.indexOf(["org.clothocad.model.NucSeq","Part","Vector"],j(a))>0:j(a)==b}function s(a){return M[a]||M["default"]}function t(a){var b;return b=B[a]?angular.extend({schema:B[a].schema},B[a].scaffold):angular.extend({schema:a},B.Instance.scaffold)}function u(a){return angular.map(C,function(b,c){return a[c]})}var v=new b("clothoSchemas","#992299"),w="org.clothocad.core.schema.Schema",x="org.clothocad.core.schema.BuiltInSchema",y="org.clothocad.core.schema.ClothoSchema",z=d.defer(),A={},B={Instance:{readable:"Instance",editor_template_url:"views/_editor/sharable.html",schema:!1,scaffold:{language:"JSONSCHEMA"},"class":"info"},Function:{readable:"Function",editor_template_url:"views/_editor/function.html",schema:"org.clothocad.core.datums.Function",scaffold:{},"class":"success"},Schema:{readable:"Schema",editor_template_url:"views/_editor/schema.html",schema:y,scaffold:{language:"JSONSCHEMA"},"class":"danger"},View:{readable:"View",editor_template_url:"views/_editor/view.html",schema:"org.clothocad.core.datums.View",scaffold:{},"class":"warning"}},C={name:{name:"name",type:"string",description:"Name of the Instance, given by author"},id:{name:"id",type:"string",description:"Unique ID referring to this object"},schema:{name:"schema",type:"string",description:"Pattern describing contents and organization of instance data"},description:{name:"description",type:"string",description:"Description of object, written by author"},author:{name:"author",type:"string",description:"User who created this object"}},D=[{name:"Public",value:"PUBLIC"},{name:"Private",value:"PRIVATE"},{name:"Read Only",value:"READONLY"}],E=[{name:"RegExp",value:"regex"},{name:"Not Null",value:"notnull"}],F={"boolean":"java.lang.Boolean",number:"java.lang.Long",object:"java.util.HashMap",array:"java.util.List",string:"java.lang.String"},G=/[^\s"']+|"([^"]*)"|'([^']*)'/,H=(G.toString(),{"boolean":{type:"boolean",input:{type:"checkbox"}},number:{type:"number",input:{type:"number"}},string:{type:"string",input:{type:"text"}},date:{type:"date",input:{type:"date"}},array:{type:"array"},"org.bson.types.ObjectId":{type:"id",input:{type:"text"}},object:{type:"object"},"default":{type:"string",input:{type:"text"}}}),I={"java.lang.Boolean":{type:"boolean",input:{type:"checkbox"}},"java.lang.Number":{type:"number",input:{type:"number"}},"java.lang.Long":{type:"number",input:{type:"number"}},"java.lang.Integer":{type:"number",input:{type:"number","ng-pattern":"/^[-+]?\\d+$/"}},"java.lang.Short":{type:"number",input:{type:"number",min:-32767,max:32767}},"java.lang.String":{type:"string",input:{type:"text"}},"java.awt.Color":{type:"color",input:{type:"color"}},"java.util.Date":{type:"date",input:{type:"date"}},"java.util.Set":{type:"array"},"java.util.List":{type:"array"},"org.bson.types.ObjectId":{type:"id",input:{type:"text"}},"java.util.Map":{type:"object"},"java.util.HashMap":{type:"object"},"default":{type:"string",input:{type:"text"}}};a.query({schema:w},{mute:!0}).then(function(a){angular.forEach(a,function(a){A[a.id]=a.name}),z.resolve(a)});var J=function(b){return a.run("clotho.functions.schema.getParents",[b],{mute:!0})},K=function(a){function b(a){a.superClass?e.then(function(){return h(a).then(function(a){return c=c.concat(a.fields),b(a)},function(){v.warn("couldnt get parent schema "+a.superClass),f.reject(c)})}):f.resolve()}var c=[];if(angular.isUndefined(a)||!a.superClass)return d.when(c);var e=d.when(),f=d.defer();return b(a),f.promise.then(function(){return e}).then(function(){return c})},L=function(b){function c(b){b.superClass?f.then(function(){return a.get(b.superClass).then(function(a){return e.fields=e.fields.concat(a.fields),c(a)},function(){v.warn("couldnt get parent schema "+b.superClass),g.reject(e)})}):g.resolve()}if(angular.isUndefined(b))return d.when();if(!b.superClass)return d.when(b);var e=angular.copy(b),f=d.when(),g=d.defer();return c(b),g.promise.then(function(){return f}).then(function(){return e})},M={Instance:"glyphicon glyphicon-file",Function:"glyphicon glyphicon-play-circle",View:"glyphicon glyphicon-picture",Schema:"glyphicon glyphicon-cog","default":"glyphicon glyphicon-file"};return{retrievedSchemas:z.promise,downloadSchemaDependencies:L,getSuperclassFields:K,getParentSchemaIds:J,sharableTypes:B,isBasicField:function(a){return angular.isDefined(C[angular.lowercase(a)])},isPrimitiveField:function(a){return angular.isDefined(F[angular.lowercase(a)])},accessTypes:D,constraintTypes:E,primitiveToJava:F,formTypeMap:H,javaToJavascript:I,isSharable:i,isFunction:n,isView:o,isSchema:k,isInstance:function(a){return"Instance"==q(a)},isBuiltIn:l,isClothoSchema:m,mapSchemaIdToName:g,isInstanceOfSchema:r,dirtyDetermineType:q,determineSharableType:p,determineSharableIcon:s,determineFieldType:e,determineSchema:j,createScaffold:t,typeToColorClass:f,sharableBasicFields:C,pruneToBasicFields:u}}]),angular.module("clotho.clothoDirectives").directive("clothoRun",["Clotho","$parse",function(a,b){var c={input:!0,textarea:!0,select:!0};return{restrict:"A",require:"ngModel",scope:!0,link:function(d,e,f,g){function h(a){angular.isDefined(l.assign)&&l.assign(d,a)}function i(a){var b=k?"val":"text";e[b](a)}function j(a){return angular.isArray(a)?a:[a]}var k=!!c[angular.lowercase(e[0].nodeName)];k&&(g.$render=angular.noop);var l=b(f.clothoRunUpdate);d.$watch(function(){return f.clothoRun},function(a){a&&m(g.$modelValue)}),d.$watch(function(){return g.$modelValue},function(a){m(a)});var m=function(b){return b=j(b),a.run(f.clothoRun,b).then(function(a){h(a),i(a)})}}}}]),angular.module("clotho.foundation").directive("clothoTool",["$http","$compile","Debug",function(a,b,c){var d=new c("clothoTool","#cc9999"),e="partials/tools/",f={revcomp:{partial:"revcomp.html"},digestCuts:{partial:"digestCuts.html",dependencies:{mixin:e+"digestCuts.js"}},ligation:{partial:"ligation.html",dependencies:{mixin:e+"ligation.js"}},pcr:{partial:"pcr.html",dependencies:{mixin:e+"pcr.js"}}};return{restrict:"A",replace:!0,link:function(c,g,h){function i(i){$clotho.extensions.downloadDependencies(f[i].dependencies).then(function(j){a.get(e+f[i].partial).success(function(a){g.removeClass("loading"),g.html(b(a)(c)),j()}).error(function(){d.error("Could not find tool"+h.clothoTool),g.remove()})})}return angular.isUndefined(h.clothoTool)?void d.warn("clothoTool attr is not defined"):(g.addClass("clotho-tool loading"),void c.$watch(function(){return h.clothoTool},function(a){f[a]?i(a):d.warn("tool name not provided")}))}}}]),angular.module("clotho.clothoDirectives").directive("clothoShow",["$q","$http","$timeout","$browser","$rootScope","$compile","Clotho","PubSub","ClothoUtils",function(a,b,c,d,e,f,g,h,i){var j=i.generateWidgetUrl,k=function(a){return b.get(j(a)+"/model.json").then(function(a){return a.data})};return{terminal:!0,restrict:"E",scope:{id:"@",callback:"=?"},controller:["$scope","$element","$attrs",function(){}],link:function(b,e,g){console.log("directive linked"),b.$watch("id",function(l){l&&(console.log(l),e.addClass("clothoWidget"),a.when(k(b.id)).then(function(a){return i.downloadViewDependencies(a)}).then(function(a){if(a.dictionary=angular.extend({},a.dictionary,a.importedViews),a.dictionary.id=a.id,a.bootstrap){var i=a.id+"-additions",k=a.bootstrap.excludeExtensionsModule!==!1?["clotho.extensions"]:[];angular.module(i,k).run(["$rootScope",function(b){angular.extend(b,a.dictionary),b.prefixUrl=function(b,c){return j(c?c:a.id,b)}}]);var l=[];l.push(function(a,b){a.value("$anchorScroll",angular.noop),a.value("$browser",d),a.service("lazyScripts",["$q","$timeout","$document",function(a,b,c){var d=[];this.register=function(a){d.push($clotho.extensions.mixin(a))},b(function(){a.all(d).then(function(){c.triggerHandler("WidgetContentLoaded")})})}]),b.directive("script",["$parse","$rootScope","lazyScripts",function(a,b,c){return{restrict:"E",terminal:!0,compile:function(d,e){if(e.ngSrc){var f=a(e.ngSrc)(b);c.register(f)}}}}])}),l=l.concat(a.bootstrap.modules,i),e.html("<div ng-include=\"prefixUrl('index.html')\"></div>"),a.controller&&g.$set("ng-controller",a.controller),e.data("$injector",null),angular.bootstrap(e,l)}else{var m;if(_.indexOf(a.files,"index.html")>=0)m="<div ",a.controller&&(m+='ng-controller="'+a.controller+'" '),m+="ng-include=\"prefixUrl('index.html')\"></div>";else if(1==a.files.length){var n=/([a-z\-_0-9\/\:\.]*\.(jpg|jpeg|png|gif))/i;n.test(a.files[0])&&(m='<img src="'+j(a.id,a.files[0])+'"alt="view '+a.id+'" />')}angular.extend(b,a.dictionary),e.html(f(m)(b))}c(function(){angular.isFunction(b.callback)&&b.callback(e),h.trigger("clothoShow:"+b.id,[b.id,e,a])})}))})}}}]),angular.module("clotho.clothoDirectives").directive("clothoModal",["$parse","$timeout","$http","$compile","$sce","$injector",function(a,b,c,d,e,f){var g,h=f.has("hotkeys");return h&&(g=f.get("hotkeys")),{restrict:"E",replace:!0,transclude:!0,templateUrl:"views/_foundation/clothoModal.html",scope:{id:"@?",open:"=?",onClose:"&?",onOpen:"&?",title:"@?",content:"=?",templateUrl:"=?",actions:"=?",modalSize:"@"},controller:["$scope","$element","$attrs",function(a){a.$close=function(){a.open&&(a.open=!1,g.del("esc"),b(function(){a.onClose()}))}}],link:function(a,f,i){i.open||b(function(){a.open=!0}),a.$watch("content",function(b){a.contentTrusted=e.trustAsHtml(b)}),a.$watch("templateUrl",function(b,e){!b||b==e&&e&&a.hasTemplate||c.get(b,{cache:!0}).success(function(b){angular.element(f[0].querySelector("[template-insert]")).html(d(b)(a)),a.hasTemplate=!0}).error(function(){a.hasTemplate=!1})}),a.$watch("open",function(b){b&&(a.open=!0,h&&g.bindTo(a).add("esc",a.$close),angular.isFunction(a.onOpen)&&a.onOpen())}),a.$on("$destroy",function(){b(function(){a.onClose()})})}}}]).service("$clothoModal",["$window","$rootScope","$compile","$parse",function(a,b,c,d){function e(){i(h),h&&h.$destroy(),g&&g.remove(),g=null,h=null,i=angular.noop}var f=angular.element(a.document.body),g=null,h=null,i=angular.noop;this.create=function(a,i){e();var j=angular.extend({},a);h=angular.isDefined(i)?i.$new():b.$new();var k=d(j["on-close"])||angular.noop;h.clothoModalClose=function(){k(h),e()},j["on-close"]="clothoModalClose()",g=c(angular.element("<clotho-modal>").attr(j))(h),f.append(g)},this.destroy=e}]),angular.module("clotho.clothoDirectives").provider("$clothoPopup",function(){function a(a){var b=a||"click",d=c[b]||b;return{show:b,hide:d}}function b(a){var b=/[A-Z]/g,c="-";return a.replace(b,function(a,b){return(b?c:"")+a.toLowerCase()})}var c={mouseenter:"mouseleave",click:"click",focus:"blur"},d={};this.$get=["$animate","$injector","$window","$document","$compile","$timeout","$parse","Clotho","ClothoSchemas",function(c,e,f,g,h,i){var j,k=e.has("hotkeys");return k&&(j=e.get("hotkeys")),function(e,l){l=angular.extend({},d,l);var m=g.find("body"),n=function(a){var b=a[0].getBoundingClientRect();return{width:b.width||a.prop("offsetWidth"),height:b.height||a.prop("offsetHeight"),top:b.top+(f.pageYOffset||g[0].documentElement.scrollTop),left:b.left+(f.pageXOffset||g[0].documentElement.scrollLeft)}},o="<div "+b(e)+'-inner sharable-id="sharable_id" sharable-model="passedModel" popup-title="popup_title" popup-placement="{{popup_placement}}" reposition="repositionFunction()"></div>';return{restrict:"EA",scope:{passedModel:"=?"+e+"Model"},compile:function(){var b=h(o);return function(d,f,g){function h(){var a,b,c,e;switch(a=n(f),b=x.prop("offsetWidth"),c=x.prop("offsetHeight"),d.popup_placement){case"topRight":e={top:a.top+a.height/2-55,left:a.left+a.width};break;case"topLeft":e={top:a.top+a.height/2-55,left:a.left-b};break;case"right":e={top:a.top+a.height/2-c/2,left:a.left+a.width};break;case"top":e={top:a.top-c,left:a.left+a.width/2-b/2};break;case"left":e={top:a.top+a.height/2-c/2,left:a.left-b};break;case"bottom":e={top:a.top+a.height,left:a.left+a.width/2-b/2};break;default:e={top:a.top+a.height,left:a.left}}e.top+="px",e.left+="px",x.toggleClass("right","topRight"==d.popup_placement),x.toggleClass("left","topLeft"==d.popup_placement),x.css(e)}function o(){return d.popupForceOpen||d.popupOpen}function p(){d.popupOpen?r():q()}function q(){t()()}function r(){d.$apply(function(){u()})}function s(a){angular.isDefined(a)?a?t()():u():o()?t()():u()}function t(){return angular.isEmpty(d.sharable_id)&&angular.isEmpty(d.passedModel)?(d.popupOpen=!1,angular.noop):(v(),x.css({top:0,left:0,display:"block"}),c.enter(x,m,angular.element(m[0].lastChild),angular.noop),h(),d.popupOpen=!0,k&&j.bindTo(y).add("esc",u),h)}function u(a){(a||!d.popupForceOpen)&&(d.popupOpen=!1,x&&w())}function v(){x&&w(),y=d.$new(),x=b(y,function(){}),y.$digest()}function w(){x&&(x.remove(),x=null),y&&(y.$destroy(),y=null)}if(!angular.isUndefined(g[e+"Id"])||!angular.isUndefined(g[e+"Model"])){var x,y,z,A=!1;d.repositionFunction=function(){i(function(){h()})},g.$observe(e+"Id",function(a,b){!a||b&&a==b||(d.sharable_id=a)}),g.$observe(e+"Title",function(a,b){!a||b&&a==b||(d.popup_title=a)}),g.$observe(e+"Placement",function(a){d.popup_placement=a||l.placement});var B=function(){A&&(f.unbind(z.show,q),f.unbind(z.hide,r))};g.$observe(e+"Trigger",function(b){B(),"none"!=b&&(z=a(b||l.trigger),z.show===z.hide?f.bind(z.show,p):(f.bind(z.show,q),f.bind(z.hide,r)),A=!0)}),d.$watch(function(){return g[e+"Open"]},function(a){d.popupForceOpen=d.$eval(a),d.popupOpen=d.$eval(a),setTimeout(function(){s(d.popupOpen)})}),d.$on("$locationChangeSuccess",function(){d.popupOpen&&u(!0)}),d.$on("$destroy",function(){B(),u(!0)})}}}}}}]}),angular.module("clotho.clothoDirectives").directive("sharablePopup",["$clothoPopup",function(a){return a("sharablePopup")}]).directive("sharablePopupInner",["Clotho","ClothoSchemas","$timeout","$injector",function(a,b,c,d){var e=d.has("clothoEditorDirective");return{restrict:"EA",replace:!0,scope:{sharableId:"=?",sharableModel:"=?",placement:"@popupPlacement",reposition:"&"},templateUrl:"views/_foundation/sharableBasicFieldsPopup.html",link:function(d){function f(c){d.sharable=c,b.determineSharableType(c).then(function(a){d.type=a},function(){d.type=b.dirtyDetermineType(c)}).then(function(){d.iconClass=b.determineSharableIcon(d.type),d.labelClass="label-"+b.typeToColorClass(d.type)}),b.isSchema(c)&&(d.isSchema=!0,a.get(c.id,{mute:!0}).then(function(a){b.getSuperclassFields(a).then(function(b){d.schema=a,d.inheritedFields=b})}))}d.editorPresent=e,d.$watch("sharableModel",function(b){b&&(f(b),b.id&&a.get(b.id,{mute:!0}).then(function(a){d.fullSharable=a}))}),d.$watch("sharableId",function(c){c&&a.get(c,{mute:!0}).then(function(a){angular.isEmpty(a)||(d.fullSharable=a,f(b.pruneToBasicFields(a)),d.reposition())})}),d.showView=function(){d.activeView=d.activeView?"":d.type,c(function(){d.reposition()})},d.edit=a.edit,d.$on("$destroy",function(){})}}}]),!function(a,b,c){function d(a,c){var d=b.createElement("script"),e=k;d.onload=d.onerror=d[p]=function(){d[n]&&!/^c|loade/.test(d[n])||e||(d.onload=d[p]=null,e=1,c())},d.async=1,d.src=a,f.insertBefore(d,f.firstChild)}function e(a,b){q(a,function(a){return!b(a)})}var f=b.getElementsByTagName("head")[0],g={},h={},i={},j={},k=!1,l="push",m="DOMContentLoaded",n="readyState",o="addEventListener",p="onreadystatechange",q=function(a,b){for(var c=0,d=a.length;d>c;++c)if(!b(a[c]))return k;return 1};!b[n]&&b[o]&&(b[o](m,function t(){b.removeEventListener(m,t,k),b[n]="complete"},k),b[n]="loading");var r=function(a,b,f){function k(){if(!--s){g[p]=1,o&&o();for(var a in i)q(a.split("|"),m)&&!e(i[a],m)&&(i[a]=[])}}function m(a){return a.call?a():g[a]}a=a[l]?a:[a];var n=b&&b.call,o=n?b:f,p=n?a.join(""):b,s=a.length;return c(function(){e(a,function(a){j[a]?(p&&(h[p]=1),k()):(j[a]=1,p&&(h[p]=1),d(r.path?r.path+a+".js":a,k))})},0),r};r.get=d,r.ready=function(a,b,c){a=a[l]?a:[a];var d=[];return!e(a,function(a){g[a]||d[l](a)})&&q(a,function(a){return g[a]})?b():!function(a){i[a]=i[a]||[],i[a][l](b),c&&c(d)}(a.join("|")),r};var s=a.$script;r.noConflict=function(){return a.$script=s,this},"undefined"!=typeof module&&module.exports?module.exports=r:a.$script=r}(this,document,setTimeout),angular.module("clotho.extensions",[]).config(["$controllerProvider","$compileProvider","$filterProvider","$provide",function(a,b,c,d){window.$clotho.extensions=$clotho.extensions={},$clotho.extensions.providers={$controllerProvider:a,$compileProvider:b,$filterProvider:c,$provide:d},$clotho.extensions._controller=$clotho.extensions.controller,$clotho.extensions._service=$clotho.extensions.service,$clotho.extensions._factory=$clotho.extensions.factory,$clotho.extensions._value=$clotho.extensions.value,$clotho.extensions._directive=$clotho.extensions.directive,$clotho.extensions.controller=function(b,c){return a.register(b,c),this},$clotho.extensions.service=function(a,b){return d.service(a,b),this},$clotho.extensions.factory=function(a,b){return d.factory(a,b),this},$clotho.extensions.value=function(a,b){return d.value(a,b),this},$clotho.extensions.directive=function(a,c){return b.directive(a,c),this},$clotho.extensions.filter=function(a,b){return c.filter(a,b),this}}]).run(["$rootScope","$q","$timeout","$templateCache","$http","$rootElement","$compile",function(a,b,c,d,e,f){{var g=function(){return angular.module("clotho.extensions")._invokeQueue};g().length}$clotho.extensions.recompile=function(a,b){"undefined"!=typeof a&&(b=b||{},a.hasClass("ng-scope")&&a.scope().$destroy(),f.injector().invoke(function(c,d){var e=d.$new(!0);angular.extend(e,b),c(a)(e),d.$apply()}))},$clotho.extensions.extendPrimaryRootscope=function(b){$clotho.extensions.extend(a,b)},$clotho.extensions.downloadDependencies=function(a){return angular.isEmpty(a)?b.when(angular.noop):b.all([$clotho.extensions.css(a.css),$clotho.extensions.mixin(a.mixin),$clotho.extensions.script(a.script)]).then(function(){return function(){c(function(){$clotho.extensions.script(a.onload)})}})},$clotho.extensions.mixin=function(d){if(angular.isUndefined(d)||""==d)return b.when("no mixin url");var e=b.defer(),f=c(function(){e.reject(null)},5e3);return $script(d,function(){c.cancel(f),a.$safeApply(e.resolve(d))}),e.promise},$clotho.extensions.script=function(a){if(angular.isUndefined(a)||0==a.length)return b.when("no script url");var c=[];return angular.isString(a)&&(a=[a]),angular.forEach(a,function(a){c.push(a+"?_="+Date.now())}),$clotho.extensions.mixin(c)};var h=[];$clotho.extensions.css=function(d){if(angular.isUndefined(d)||""==d)return b.when("no css url");if(_.indexOf(h,d)>-1)return b.when("CSS url already added");var e=b.defer(),f=c(function(){e.reject(null)},5e3);if(document.createStyleSheet)document.createStyleSheet(d),a.$safeApply(e.resolve());else{var g=document.createElement("link");g.type="text/css",g.rel="stylesheet",g.href=d,document.getElementsByTagName("head")[0].appendChild(g),c.cancel(f),a.$safeApply(e.resolve())}return h.push(d),e.promise},$clotho.extensions.cache=function(a){if(angular.isUndefined(a)||""==a)return b.when();var f=b.defer(),g=c(function(){f.reject(null)},5e3);return e.get(a).success(function(b){c.cancel(g),d.put(a,b),f.resolve(b)}).error(function(a){f.reject(a)}),f.promise},$clotho.extensions.bootstrap=angular.bootstrap,$clotho.extensions.determineUrlExtension=function(a){var b=a.split("?")[0];return b.substr(b.lastIndexOf(".")+1)};document.getElementsByTagName("head")[0]}]),angular.module("clotho.foundation").directive("clothoShowAuth",["PubSub",function(a){function b(a,b){var c=a.$eval(b);return"string"==typeof c||angular.isArray(c)||(c=b),"string"==typeof c&&(c=c.split(",")),c}function c(a,b){var c=!1;return angular.forEach(b,function(b){return b===a?(c=!0,!0):!1}),c}function d(a){if(!a.length)throw new Error("clotho-show-auth directive must be login, logout, or error (you may use a comma-separated list)");return angular.forEach(a,function(a){if(!c(a,["login","logout","error"]))throw new Error('Invalid state "'+a+'" for clotho-show-auth directive, must be one of login, logout, or error')}),!0}var e="logout";return a.on("auth:login",function(){e="login"}),a.on("auth:logout",function(){e="logout"}),a.on("auth:error",function(){e="error"}),{restrict:"A",link:function(f,g,h){function i(){var a=c(e,j);setTimeout(function(){g.toggleClass("ng-cloak",!a)},0)}var j=b(f,h.clothoShowAuth);d(j),i(),a.on("auth:login auth:logout auth:error",i)}}}]),angular.module("clotho.foundation").service("Facebook",["$q","$window","$filter",function(a,b,c){function d(){g={},h=a.defer(),i=h.promise}function e(){l=!0;var b=a.defer();return o.then(function(){FB.api("/me",function(a){f(a),b.resolve(a)})}),b.promise}function f(a){l=!0,m=!1,angular.extend(g,a,{provider:"facebook",icon:"//graph.facebook.com/"+a.id+"/picture?width=500&height=500"}),h.resolve(g)}var g,h,i,j=this,k=0x52d209de19bff,l=!1,m=!1,n=a.defer(),o=n.promise;d(),b.fbAsyncInit=function(){FB.init({appId:k,xfbml:!0,version:"v2.0",cookie:!0,status:!0}),FB.getLoginStatus(function(a){"connected"===a.status?e():"not_authorized"===a.status,n.resolve(FB)}),FB.Event.subscribe("auth.login",function(){}),FB.Event.subscribe("auth.logout",function(){})},j.login=function(b){var c=angular.extend({scope:"public_profile email"},b),d=a.defer();return o.then(function(a){l?i.then(function(a){d.resolve(a)}):a.login(function(a){console.log("FB login response",a),a.authResponse?e().then(function(){i.then(function(a){d.resolve(a)})}):d.reject(a)},c)}),d.promise},j.logout=function(a){return o.then(function(b){b.getLoginStatus(function(c){"connected"===c.status&&b.logout(a)}),d()})},j.getUser=function(){var b=a.defer();return angular.isEmpty(g)?o.then(function(){l?i.then(function(a){b.resolve(a)}):b.reject()}):b.resolve(g),b.promise},j.convertToPersonSharable=function(a){return a=a||g,angular.isEmpty(g)?void 0:{schema:"org.clothocad.model.Person",id:a.email,name:a.name,emailAddress:a.email,dateCreated:c("date")(Date.now().valueOf(),"yyyy-MM-dd"),icon:a.icon,social:{facebook:a.link}}},function(a,b,c){var d,e=a.getElementsByTagName(b)[0];a.getElementById(c)||(d=a.createElement(b),d.id=c,d.src="//connect.facebook.net/en_US/sdk.js",e.parentNode.insertBefore(d,e))}(document,"script","facebook-jssdk")}]);
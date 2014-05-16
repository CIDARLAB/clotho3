angular.module("clotho.tokenizer",[]),angular.module("clotho.commandbar",["clotho.core","clotho.tokenizer","ui.keypress"]),angular.module("ui.keypress",[]).factory("keypressHelper",["$parse",function(a){var b={8:"backspace",9:"tab",13:"enter",27:"esc",32:"space",33:"pageup",34:"pagedown",35:"end",36:"home",37:"left",38:"up",39:"right",40:"down",45:"insert",46:"delete"},c=function(a){return a.charAt(0).toUpperCase()+a.slice(1)};return function(d,e,f,g){var h,i=[];h=e.$eval(g["ui"+c(d)]),angular.forEach(h,function(b,c){var d,e;e=a(b),angular.forEach(c.split(" "),function(a){d={expression:e,keys:{}},angular.forEach(a.split("-"),function(a){d.keys[a]=!0}),i.push(d)})}),f.bind(d,function(a){var c=!(!a.metaKey||a.ctrlKey),f=!!a.altKey,g=!!a.ctrlKey,h=!!a.shiftKey,j=a.keyCode;"keypress"===d&&!h&&j>=97&&122>=j&&(j-=32),angular.forEach(i,function(d){var i=d.keys[b[j]]||d.keys[j.toString()],k=!!d.keys.meta,l=!!d.keys.alt,m=!!d.keys.ctrl,n=!!d.keys.shift;i&&k===c&&l===f&&m===g&&n===h&&e.$apply(function(){d.expression(e,{$event:a})})})})}}]),angular.module("ui.keypress").directive("uiKeydown",["keypressHelper",function(a){return{link:function(b,c,d){a("keydown",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeypress",["keypressHelper",function(a){return{link:function(b,c,d){a("keypress",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeyup",["keypressHelper",function(a){return{link:function(b,c,d){a("keyup",b,c,d)}}}]),angular.module("clotho.commandbar").service("CommandBar",["Clotho","ClientAPI","Debug","$timeout","$q","$document",function(a,b,c,d,e,f){function g(a){m.unread=m.unread&&!n.log?m.unread+1:1,m.entries.unshift(a),i.log("LOG - entries: ",m.entries),n.toggle("logSnippet",!0),m.startLogTimeout()}var h={dateFilter:"mediumTime",timeFilter:"timestamp"},i=new c("Command Bar","#ffbb55"),j=function(){return angular.element(f[0].querySelector("[clotho-command-bar]"))},k=function(){return angular.element(f[0].querySelector("[clotho-command-bar] [clotho-autocomplete]"))},l=function(){k().focus()},m={};m.entries=[{text:"Welcome to Clotho!",from:"server","class":"success",timestamp:Date.now()}];var n={};n.query="",n.log=!1,n.logSnippet=!1,n.toggle=function(a,b){n[a]=angular.isDefined(b)?b:!n[a]},n.toggleActivityLog=function(){n.log=!n.log,n.log&&(m.unread="")},m.timeout=null,m.startLogTimeout=function(){m.cancelLogTimeout(),m.timeout=d(function(){n.toggle("logSnippet",!1)},1e4)},m.cancelLogTimeout=function(){d.cancel(m.timeout)};var o=function(c){if((angular.isEmpty(c)||!angular.isObject(c))&&(c=n.query),c.query=c.query.trim(),c.query){var d={"class":"info",from:"client",text:c.query,timestamp:Date.now()};return b.say(d),a.submit(c).then(function(a){console.log("resetting query"),n.query="",b.say({text:a,"class":"success"})},function(){})}return e.when(!1)};return a.listen("activityLog",function(a){g(a)},"searchbar"),{options:h,display:n,log:m,setQuery:function(a,b){angular.isDefined(b)&&b.preventDefault(),a&&(n.query=angular.isEmpty(a.value)?a.text:a.value)},submit:o,getCommandBarElement:j,getCommandBarInput:k,focusInput:l}}]),angular.module("clotho.commandbar").directive("clothoCommandBar",["Clotho","CommandBar","$location","$window","$compile","$clothoModal",function(a,b,c,d,e,f){return{restrict:"A",replace:!0,scope:{},templateUrl:"views/_command/commandbar.html",controller:["$scope","$element","$attrs",function(a){a.options=b.options,a.log=b.log,a.autocomplete=b.autocomplete,a.display=b.display,a.setQuery=b.setQuery,a.submit=b.submit,a.execute=b.execute;var c=!1;a.toggleLogin=function(a){c=angular.isDefined(a)?a:!c,c?f.create({title:"Clotho Login","template-url":"'views/_command/simpleLogin.html'"}):f.destroy()}}],link:function(b){b.showMeHow=function(){a.query({name:"Learning Clotho"}).then(function(a){c.path("/trails/"+a[0].id)})},b.goHome=function(){c.path("/")},b.aboutClotho=function(){c.path("/about")},b.teamClotho=function(){c.path("/team")}}}}]),angular.module("clotho.commandbar").controller("TerminalCtrl",["$scope","CommandBar",function(a,b){a.log=b.log}]),angular.module("clotho.commandbar").directive("logEntries",function(){return{restrict:"A",templateUrl:"views/_command/logEntries.html",scope:{entries:"=logEntries"}}}),angular.module("clotho.commandbar").controller("loginCtrl",["$scope","$timeout","Clotho",function(a,b,c){a.notification={},a.cred={username:"",password:""},a.login=function(){c.login(a.cred.username,a.cred.password).then(function(b){console.log("run login",b),b?a.notification={"class":"alert-success",message:"Log in Success"}:(a.notification={"class":"alert-danger",message:"Log in Error"},a.cred={username:"",password:""})})}}]),angular.module("clotho.tokenizer").directive("clothoAutocomplete",["Clotho","$q","$parse","$timeout","$compile","$filter","$document",function(a,b,c,d,e,f,g){var h=[8,9,13,27,37,38,39,40],i=32,j=i,k=" ",l=0;return{restrict:"A",require:"?ngModel",link:function(b,i,m){function n(){b.query=""}function o(){b.autocompletions=[],b.activeIdx=-1}function p(){n(),o(),b.hasFocus=!1,b.tokenCollection.unsetActive(),b.$digest()}function q(a){b.hasFocus&&!b.tokenCollection.isActive()&&b.activeIdx<0&&(i[0].contains(a.target)||d(p))}var r=/^['"].*/,s=c(m.autocompleteOnSelect),t=angular.element("<clotho-autocomplete-listing></clotho-autocomplete-listing>");t.attr({autocompletions:"autocompletions",active:"activeIdx",select:"select(activeIdx)","has-focus":"hasFocus",query:"query"}),b.hasFocus=!1;var u=function(c){if(r.test(c)&&(c=c.substring(1)),0!==c.length){a.autocomplete(c).then(function(a){a&&a.length&&b.query.length?b.autocompletions=f("limitTo")(a,10):o()})}};b.query="";var v;b.$watch("query",function(a){a&&a.length?(b.hasFocus=!0,b.tokenCollection.unsetActive(),l>0?(v&&d.cancel(v),v=d(function(){u(a)},l)):u(a)):o()}),b.select=function(a){var c=a>-1?b.autocompletions[a]:b.query;c&&s(b,{$item:c,$query:b.query}),o(),n(),d(function(){i[0].focus()},0,!1)},i.bind("keydown",function(a){if(a.which===j&&!r.test(b.query.charAt(0)))return void b.$apply(function(){b.select(1==b.autocompletions.length?0:-1)});if(-1!==h.indexOf(a.which)){if(8===a.which){if(b.query.length)return;if(b.tokenCollection)if(b.tokenCollection.isActive()){var c=b.tokenCollection.whichActive();b.tokenCollection.removeActiveToken(),b.tokenCollection.setActive(c)}else b.tokenCollection.setLastActive();b.$digest()}else if(40===a.which)b.autocompletions.length&&(b.activeIdx=(b.activeIdx+1)%b.autocompletions.length,b.$digest());else if(38===a.which)b.autocompletions.length&&(b.activeIdx=(b.activeIdx?b.activeIdx:b.autocompletions.length)-1,b.$digest());else if(37===a.which){if(!b.tokenCollection.isActive())return;b.tokenCollection.setPrevActive(),b.$digest()}else if(39===a.which){if(!b.tokenCollection.isActive())return;b.tokenCollection.isLastActive()?b.tokenCollection.unsetActive():b.tokenCollection.setNextActive(),b.$digest()}else 13===a.which||9===a.which?b.activeIdx>=0?b.$apply(function(){b.select(b.activeIdx)}):(b.query.length&&b.$apply(function(){b.select()}),b.$apply(function(){b.submit()})):27===a.which&&(o(),b.tokenCollection.unsetActive(),b.$digest(),b.query.length?a.stopPropagation():(b.hasFocus=!1,i.blur(),b.$digest()));a.preventDefault()}}),i.on("focus",function(){d(function(){b.hasFocus=!0})}),i.on("paste",function(){d(function(){angular.forEach(b.query.split(k),function(a){b.addToken(a)}),n(),o()})}),b.$on("$locationChangeSuccess",function(){setTimeout(p)}),g.bind("click",q),b.$on("$destroy",function(){g.unbind("click",q)}),o(),i.after(e(t)(b))}}}]),angular.module("clotho.tokenizer").directive("clothoAutocompleteListing",function(){return{restrict:"EA",scope:{autocompletions:"=",query:"=",active:"=",hasFocus:"=",select:"&"},replace:!0,templateUrl:"views/_command/autocompleteListing.html",link:function(a){a.isOpen=function(){return a.hasFocus&&a.matches.length>0},a.isActive=function(b){return a.active==b},a.selectActive=function(b){a.active=b},a.selectMatch=function(b){a.select({activeIdx:b})}}}}),angular.module("clotho.tokenizer").directive("clothoAutocompleteMatch",["ClothoSchemas",function(a){return{restrict:"EA",replace:!0,scope:{index:"=",match:"=",query:"="},templateUrl:"views/_command/autocompleteMatch.html",link:function(b,c,d){b.$watch(function(){return d.active},function(a){b.active=b.$eval(a)}),b.$watch("match",function(c){b.iconClass=a.determineSharableIcon(a.determineSharableType(c))})}}}]).filter("clothoAutocompleteHighlight",function(){function a(a){return a.replace(/([.?*+^$[\]\\(){}|-])/g,"\\$1")}return function(b,c){return c?b.replace(new RegExp(a(c),"gi"),"<strong>$&</strong>"):b}}),angular.module("clotho.tokenizer").directive("clothoToken",["Clotho","clothoTokenFactory",function(){return{restrict:"E",replace:!0,templateUrl:"views/_command/token.html",scope:{tokenCollection:"=",tokenIndex:"=",tokenActive:"=",token:"=ngModel",onRemove:"&?"},link:function(a,b){b.on("click",function(){a.tokenCollection[a.tokenActive?"unsetActive":"setActive"](a.tokenIndex)}),a.removeToken=function(b){b.preventDefault(),a.onRemove({$token:a.token,$event:b})}}}}]),angular.module("clotho.tokenizer").factory("clothoTokenCollectionFactory",["clothoTokenFactory",function(a){function b(a){this.tokens=[],this.currentTokenIndex=-1,angular.isArray(a)&&angular.forEach(a,function(a){this.addToken(a)})}return b.prototype.addToken=function(b){this.tokens.push(new a(b))},b.prototype.inRange=function(a){return a>-1&&a<this.tokens.length},b.prototype.getToken=function(a){return this.tokens[a]},b.prototype.indexOf=function(a){return this.tokens.indexOf(a)},b.prototype.removeToken=function(a){return this.inRange(a)?this.tokens.splice(a,1):!1},b.prototype.removeAll=function(){this.tokens.length=0},b.prototype.removeActiveToken=function(){if(this.isActive()){var a=this.removeToken(this.currentTokenIndex);return this.unsetActive(),a}return!1},b.prototype.setActive=function(a){return this.inRange(a)?(this.currentTokenIndex=a,a):!1},b.prototype.setLastActive=function(){this.setActive(this.tokens.length-1)},b.prototype.setPrevActive=function(){this.currentTokenIndex=(this.currentTokenIndex>0?this.currentTokenIndex:this.tokens.length)-1},b.prototype.setNextActive=function(){this.currentTokenIndex=(this.currentTokenIndex+1)%this.tokens.length},b.prototype.unsetActive=function(){this.currentTokenIndex=-1},b.prototype.isActive=function(a){return angular.isDefined(a)?this.currentTokenIndex==a:this.currentTokenIndex>-1},b.prototype.whichActive=function(){return this.currentTokenIndex},b.prototype.isLastActive=function(){return this.currentTokenIndex===this.tokens.length-1},b}]),angular.module("clotho.tokenizer").factory("clothoTokenFactory",["Clotho",function(a){function b(b){var c=this;c.model=b,this.isSharable()&&(c.fullSharablePromise=a.get(c.model.id,{mute:!0}).then(function(a){c.fullSharable=a}))}return b.prototype.readable=function(){return this.model.name||this.model},b.prototype.id=function(){return this.model.id||null},b.prototype.isAmbiguous=function(){return angular.isArray(this.model)},b.prototype.isSharable=function(){return!this.isAmbiguous()&&angular.isDefined(this.model.id)},b}]),angular.module("clotho.tokenizer").directive("clothoTokenizer",["$parse","clothoTokenCollectionFactory","Debug",function(a,b,c){var d=new c("clothoTokenizer","#ee7711");return{restrict:"E",replace:!0,require:"ngModel",templateUrl:"views/_command/tokenizer.html",link:function(c,e,f,g){function h(){d.log("updating model (query, tokens)",j,c.tokenCollection.tokens),g.$setViewValue({query:j,tokens:c.tokenCollection.tokens})}c.placeholder=f.placeholder;var i=a(f.startingTags)(c),j="";c.tokenCollection=new b(i),g.$render=function(){console.log("rendering")},c.$watchCollection("tokenCollection.tokens",function(){j="",angular.forEach(c.tokenCollection.tokens,function(a){j+=a.readable()+" "}),h()}),c.addToken=function(a){d.log("TOKENIZER_LINK adding token",a),c.tokenCollection.addToken(a)},c.removeToken=function(a){d.log("TOKENIZER_LINK removing token",a),c.tokenCollection.removeToken(a)},c.tokenActive=function(a){return c.tokenCollection.isActive(a)},c.focusInput=function(){e[0].querySelector("[clotho-autocomplete]").focus()}}}}]);
angular.module("clotho.commandbar",["clotho.core","ui.keypress"]),angular.module("ui.keypress",[]).factory("keypressHelper",["$parse",function(a){var b={8:"backspace",9:"tab",13:"enter",27:"esc",32:"space",33:"pageup",34:"pagedown",35:"end",36:"home",37:"left",38:"up",39:"right",40:"down",45:"insert",46:"delete"},c=function(a){return a.charAt(0).toUpperCase()+a.slice(1)};return function(d,e,f,g){var h,i=[];h=e.$eval(g["ui"+c(d)]),angular.forEach(h,function(b,c){var d,e;e=a(b),angular.forEach(c.split(" "),function(a){d={expression:e,keys:{}},angular.forEach(a.split("-"),function(a){d.keys[a]=!0}),i.push(d)})}),f.bind(d,function(a){var c=!(!a.metaKey||a.ctrlKey),f=!!a.altKey,g=!!a.ctrlKey,h=!!a.shiftKey,j=a.keyCode;"keypress"===d&&!h&&j>=97&&122>=j&&(j-=32),angular.forEach(i,function(d){var i=d.keys[b[j]]||d.keys[j.toString()],k=!!d.keys.meta,l=!!d.keys.alt,m=!!d.keys.ctrl,n=!!d.keys.shift;i&&k===c&&l===f&&m===g&&n===h&&e.$apply(function(){d.expression(e,{$event:a})})})})}}]),angular.module("ui.keypress").directive("uiKeydown",["keypressHelper",function(a){return{link:function(b,c,d){a("keydown",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeypress",["keypressHelper",function(a){return{link:function(b,c,d){a("keypress",b,c,d)}}}]),angular.module("ui.keypress").directive("uiKeyup",["keypressHelper",function(a){return{link:function(b,c,d){a("keyup",b,c,d)}}}]),angular.module("clotho.commandbar").service("CommandBar",["Clotho","ClientAPI","Debug","$timeout","$q","$document",function(a,b,c,d,e,f){function g(){q.log=!0}function h(a){o.unread=o.unread&&!q.log?o.unread+1:1,o.entries.unshift(a),j.log("LOG - entries: ",o.entries),q.show("logSnippet"),o.startLogTimeout()}var i={dateFilter:"mediumTime",timeFilter:"timestamp"},j=new c("Command Bar","#ffbb55"),k=function(){return angular.element(f[0].querySelector("[clotho-command-bar]"))},l=function(){return angular.element(f[0].getElementById("clotho_command_input"))},m=function(){return angular.element(f[0].getElementById("clotho_logButton"))},n=function(){l().focus()},o={},p={};p.autocompletions=[],p.autoDetail={},p.detailTemplate={},p.detailModel={},p.detailUUID=-1,o.entries=[{text:"Welcome to Clotho!",from:"server","class":"success",timestamp:Date.now()}];var q={};q.query="",q.queryHistory=[],q.autocomplete=!1,q.autocompleteDetail=!1,q.autocompleteDetailInfo=!1,q.help=!1,q.log=!1,q.logSnippet=!1,q.show=function(a){q[a]||(q[a]=!0)},q.hide=function(a){q[a]&&(q[a]=!1)},q.toggle=function(a){q[a]=!q[a]},q.genLogPos=function(){var a=m()[0];q.logpos={left:a.offsetLeft+a.scrollWidth/2-200+"px",top:a.offsetTop+a.scrollHeight+"px"}},q.detail=function(b){"undefined"!=typeof b&&(b!=p.detailUUID&&(q.hide("autocompleteDetailInfo"),p.detailUUID=b),a.autocompleteDetail(p.detailUUID).then(function(a){p.autoDetail=a,q.show("autocompleteDetail")}))},q.undetail=function(){q.hide("autocompleteDetail"),q.hide("autocompleteDetailInfo"),p.detailModel={}},q.detailInfo=function(a,b){switch(a){case"command":p.detailTemplate="views/_command/detail-command.html";break;case"author":p.detailTemplate="views/_command/detail-author.html"}p.detailModel=p.autoDetail.sharables[b],"author"==a&&(p.detailModel=p.detailModel.author),q.show("autocompleteDetailInfo")},o.timeout=null,o.startLogTimeout=function(){o.cancelLogTimeout(),o.timeout=d(function(){q.hide("logSnippet")},1e4)},o.cancelLogTimeout=function(){d.cancel(o.timeout)};var r=function(a){j.log("execute called (not implemented) "+a),q.hide("autocomplete"),q.undetail()},s=function(c){if(c||(c=q.query),c){var d={"class":"info",from:"client",text:c,timestamp:Date.now()};return q.queryHistory.push(d),q.undetail(),b.say(d),a.submit(c).then(function(a){q.query="",b.say({text:a,"class":"success"})},function(){})}return e.when(!1)};return a.listen("activityLog",function(a){h(a)},"searchbar"),{options:i,display:q,log:o,setQuery:function(a,b){"undefined"!=typeof b&&b.preventDefault(),a&&(q.query=a.value?a.value:a.text)},autocomplete:p,submit:s,execute:r,getCommandBarElement:k,getCommandBarInput:l,focusInput:n,showActivityLog:g}}]),angular.module("clotho.commandbar").directive("clothoCommandBar",["Clotho","CommandBar","$location","$window","hotkeys",function(a,b,c,d){return{restrict:"A",replace:!0,templateUrl:"views/_command/commandbar.html",controller:["$scope","$element","$attrs",function(d){d.options=b.options,d.log=b.log,d.autocomplete=b.autocomplete,d.display=b.display,d.setQuery=b.setQuery,d.submit=b.submit,d.execute=b.execute,d.$watch("display.query",function(b){d.display.autocomplete=!!b,b&&angular.isString(b)&&a.autocomplete(d.display.query).then(function(a){d.autocomplete.autocompletions=a})}),d.activeIndex=!1,d.selectAutoNext=function(a){a.preventDefault(),d.activeIndex=(d.activeIndex+1)%d.display.queryHistory.length,b.setQuery(d.display.queryHistory[d.activeIndex])},d.selectAutoPrev=function(a){a.preventDefault(),d.activeIndex=(d.activeIndex?d.activeIndex:d.display.queryHistory.length)-1,b.setQuery(d.display.queryHistory[d.activeIndex])},d.fullPageLog=function(){c.path("/terminal"),d.display.hide("log")},d.hideAutocomplete=function(){d.display.hide("autocomplete"),d.display.undetail()},d.pathIsTerminal=function(){var a=/^\/terminal.*$/;return a.test(c.path())},d.showClothoLoginModal=!1,d.showLogin=function(){d.showClothoLoginModal=!0}}],link:function(b){b.newPage=function(){d.open(d.location.origin,"_blank")},b.newWorkspace=function(){d.open(d.location.origin,"_blank")},b.showMeHow=function(){a.query({name:"Learning Clotho"}).then(function(a){c.path("/trails/"+a[0].id)})},b.goHome=function(){c.path("/")},b.aboutClotho=function(){c.path("/about")},b.teamClotho=function(){c.path("/team")},b.toggleTooltips=function(){console.log("tooltips")}}}}]),angular.module("clotho.commandbar").controller("TerminalCtrl",["$scope","CommandBar",function(a,b){a.log=b.log}]),angular.module("clotho.commandbar").directive("commandBarAutocomplete",["Clotho","CommandBar","$location","$window",function(){return{restrict:"A",templateUrl:"views/_command/autocomplete.html"}}]),angular.module("clotho.commandbar").directive("logEntries",function(){return{restrict:"A",templateUrl:"../views/_command/logEntries.html",scope:{entries:"=logEntries"}}}),angular.module("clotho.commandbar").factory("clothoTokenFactory",["Clotho",function(a){function b(b,c){this.value=b,this.uuid=c||void 0,this.isSharable=angular.isDefined(c),this.isSharable&&(this.fullSharablePromise=a.get(this.uuid).then(function(a){this.fullSharable=a}))}return b.prototype.isAmbiguous=function(){},b.prototype.isValid=function(){},b}]).factory("clothoTokenCollectionFactory",["clothoTokenFactory",function(a){function b(a){this.tokens=[],this.currentSelectedIndex=-1,angular.isArray(a)&&angular.forEach(a,function(a){this.addToken(a)})}return b.prototype.addToken=function(){this.tokens.push(new a(arguments))},b.prototype.inRange=function(a){return a>-1&&a<this.tokens.length},b.prototype.getToken=function(a){return this.tokens[a]},b.prototype.indexOf=function(a){return this.tokens.indexOf(a)},b.prototype.removeToken=function(a){return this.inRange(a)?this.tokens.splice(a,1):!1},b.prototype.removeAll=function(){this.tokens.length=0},b.prototype.removeActiveToken=function(){if(this.isActive()){var a=this.removeToken(this.currentSelectedIndex);return this.unsetActive(),a}return!1},b.prototype.setActive=function(a){return this.inRange(a)?(this.currentSelectedIndex=a,a):!1},b.prototype.setLastActive=function(){this.setActive(this.tokens.length-1)},b.prototype.setPrevActive=function(){this.currentSelectedIndex=(this.currentSelectedIndex>0?this.currentSelectedIndex:this.tokens.length)-1},b.prototype.setNextActive=function(){this.currentSelectedIndex=(this.currentSelectedIndex+1)%this.tokens.length},b.prototype.unsetActive=function(){this.currentSelectedIndex=-1},b.prototype.isActive=function(a){return angular.isDefined(a)?this.currentSelectedIndex==a:this.currentSelectedIndex>-1},b}]).directive("clothoTokenizer",["$parse","clothoTokenCollectionFactory",function(a,b){return{restrict:"E",replace:!0,require:"ngModel",templateUrl:"views/_command/tokenizer.html",controller:["$scope","$element","$attrs",function(){}],link:function(c,d,e,f){function g(){console.log("updating model",c.tokenCollection.tokens),f.$setViewValue(c.tokenCollection.tokens),console.log(f)}c.placeholder=e.placeholder;var h=a(e.startingTags)(c);c.tokenCollection=new b(h),c.$watchCollection("tokenCollection.tokens",function(){console.log("COLLECTION CHANGED"),g()}),c.addToken=function(a){console.log("TOKENIZER_LINK adding token",a),c.tokenCollection.addToken(a)},c.removeToken=function(a){console.log("TOKENIZER_LINK removing token",a),c.tokenCollection.removeToken(a)},c.tokenActive=function(a){return c.tokenCollection.isActive(a)},c.focusInput=function(){d[0].querySelector(".clothoAutocomplete").focus()}}}}]).directive("clothoAutocomplete",["Clotho","$q","$parse","$timeout","$compile","$filter",function(a,b,c,d,e,f){var g=[8,9,13,27,37,38,39,40];return{restrict:"A",controller:["$scope","$element","$attrs",function(){}],link:function(b,h,i){var j=c(i.autocompleteOnSelect),k=angular.element("<clotho-autocomplete-listing></clotho-autocomplete-listing>");k.attr({matches:"queryResults",active:"activeIdx",select:"select(activeIdx)",hasFocus:"hasFocus",query:"query"}),b.hasFocus=!1;var l=0,m=function(){b.queryResults=[],b.activeIdx=-1},n=function(c){a.autocomplete(b.query).then(function(a){a&&a.length?b.queryResults=f("limitTo")(a,10):m()})};b.query=void 0;var o;b.$watch("query",function(a){a&&a.length?(b.hasFocus=!0,b.tokenCollection.unsetActive(),l>0?(o&&d.cancel(o),o=d(function(){n(a)},l)):n(a)):m()}),b.select=function(a){var c=b.queryResults[a]||b.query;c&&j(b,{$item:c}),m(),b.query="",d(function(){h[0].focus()},0,!1)},h.bind("keydown",function(a){-1!==g.indexOf(a.which)&&(a.preventDefault(),8===a.which?b.query.length?b.$apply(function(){b.query=b.query.substring(0,b.query.length-1)}):(b.tokenCollection&&(b.tokenCollection.isActive()?b.tokenCollection.removeActiveToken():b.tokenCollection.setLastActive()),b.$digest()):40===a.which?(b.activeIdx=(b.activeIdx+1)%b.queryResults.length,b.$digest()):38===a.which?(b.activeIdx=(b.activeIdx?b.activeIdx:b.queryResults.length)-1,b.$digest()):37===a.which?(b.tokenCollection.setPrevActive(),b.$digest()):39===a.which?(b.tokenCollection.setNextActive(),b.$digest()):13===a.which||9===a.which?b.$apply(function(){b.select(b.activeIdx)}):27===a.which&&(a.stopPropagation(),m(),b.$digest()))}),h.on("blur",function(){b.hasFocus=!1,b.$apply(function(){b.tokenCollection.unsetActive()})}),h.on("focus",function(){b.hasFocus=!0}),m(),h.after(e(k)(b))}}}]).directive("clothoAutocompleteListing",function(){return{restrict:"EA",scope:{matches:"=",query:"=",active:"=",hasFocus:"=",select:"&"},replace:!0,templateUrl:"views/_command/autocompleteListing.html",link:function(a){a.isOpen=function(){return a.hasFocus&&a.matches.length>0},a.isActive=function(b){return a.active==b},a.selectActive=function(b){a.active=b},a.selectMatch=function(b){console.log(b),a.select({activeIdx:b})}}}}).directive("clothoAutocompleteMatch",function(){return{restrict:"EA",replace:!0,scope:{index:"=",match:"=",query:"="},templateUrl:"views/_command/autocompleteMatch.html",link:function(){}}}).filter("clothoAutocompleteHighlight",function(){function a(a){return a.replace(/([.?*+^$[\]\\(){}|-])/g,"\\$1")}return function(b,c){return c?b.replace(new RegExp(a(c),"gi"),"<strong>$&</strong>"):b}}).directive("clothoToken",["Clotho","clothoTokenFactory",function(){return{restrict:"E",replace:!0,templateUrl:"views/_command/token.html",scope:{tokenCollection:"=",tokenIndex:"=",tokenActive:"=",model:"=ngModel",onRemove:"&?"},controller:["$scope","$element","$attrs",function(){}],link:function(a,b){a.tokenCollection&&b.on("click",function(){a.tokenActive?(console.log("sharable object",a.fullSharable),a.tokenCollection.unsetActive(a.tokenIndex)):a.tokenCollection.setActive(a.tokenIndex)}),a.removeToken=function(b){b.preventDefault(),a.onRemove({$model:a.model})}}}}]),angular.module("clotho.commandbar").controller("loginCtrl",["$scope","$timeout","Clotho",function(a,b,c){a.notification={},a.cred={username:"",password:""},a.login=function(){c.login(a.cred.username,a.cred.password).then(function(b){console.log("run login",b),b?a.notification={"class":"alert-success",message:"Log in Success"}:(a.notification={"class":"alert-danger",message:"Log in Error"},a.cred={username:"",password:""})})}}]);
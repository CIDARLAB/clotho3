angular.module("clotho.commandbar",["clotho.core"]),angular.module("clotho.commandbar").service("CommandBar",["Clotho","ClientAPI","Debug","$timeout","$q","$document",function(a,b,c,d,e,f){function g(a){n.unread=n.unread&&!p.log?n.unread+1:1,n.entries.unshift(a),i.log("LOG - entries: ",n.entries),p.show("logSnippet"),n.startLogTimeout()}var h={dateFilter:"mediumTime",timeFilter:"timestamp"},i=new c("Command Bar","#ffbb55"),j=function(){return angular.element(f[0].getElementById("clothoCommandBar"))},k=function(){return angular.element(f[0].getElementById("clothoCommandBarInput"))},l=function(){return angular.element(f[0].getElementById("clotho_logButton"))},m=function(){j().focus()},n={},o={};o.autocompletions=[],o.autoDetail={},o.detailTemplate={},o.detailModel={},o.detailUUID=-1,n.entries=[{text:"Welcome to Clotho!",from:"server","class":"success",timestamp:Date.now()}];var p={};p.query="",p.queryHistory=[],p.autocomplete=!1,p.autocompleteDetail=!1,p.autocompleteDetailInfo=!1,p.help=!1,p.log=!1,p.logSnippet=!1,p.show=function(a){p[a]||(p[a]=!0)},p.hide=function(a){p[a]&&(p[a]=!1)},p.toggle=function(a){p[a]=!p[a]},p.genLogPos=function(){var a=l()[0];p.logpos={left:a.offsetLeft+a.scrollWidth/2-200+"px",top:a.offsetTop+a.scrollHeight+"px"}},p.detail=function(b){"undefined"!=typeof b&&(b!=o.detailUUID&&(p.hide("autocompleteDetailInfo"),o.detailUUID=b),a.autocompleteDetail(o.detailUUID).then(function(a){o.autoDetail=a,p.show("autocompleteDetail")}))},p.undetail=function(){p.hide("autocompleteDetail"),p.hide("autocompleteDetailInfo"),o.detailModel={}},p.detailInfo=function(a,b){switch(a){case"command":o.detailTemplate="views/_command/detail-command.html";break;case"author":o.detailTemplate="views/_command/detail-author.html"}o.detailModel=o.autoDetail.sharables[b],"author"==a&&(o.detailModel=o.detailModel.author),p.show("autocompleteDetailInfo")},n.timeout=null,n.startLogTimeout=function(){n.cancelLogTimeout(),n.timeout=d(function(){p.hide("logSnippet")},1e4)},n.cancelLogTimeout=function(){d.cancel(n.timeout)};var q=function(a){i.log("execute called (not implemented) "+a),p.hide("autocomplete"),p.undetail()},r=function(c){if(c||(c=p.query),c){var d={"class":"info",from:"client",text:c,timestamp:Date.now()};return p.queryHistory.push(d),p.undetail(),b.say(d),a.submit(c).then(function(a){p.query="",b.say({text:a,"class":"success"})})}return e.when(!1)};return a.listen("activityLog",function(a){g(a)},"searchbar"),{options:h,display:p,log:n,setQuery:function(a,b){"undefined"!=typeof b&&b.preventDefault(),a&&(p.query=a.value?a.value:a.text)},autocomplete:o,submit:r,execute:q,getCommandBarElement:j,getCommandBarInput:k,focusInput:m}}]),angular.module("clotho.commandbar").directive("clothoCommandBar",["Clotho","CommandBar","$location","$window",function(a,b,c,d){return{restrict:"A",replace:!0,templateUrl:"views/_command/commandbar.html",controller:["$scope","$element","$attrs",function(d){d.options=b.options,d.log=b.log,d.autocomplete=b.autocomplete,d.display=b.display,d.setQuery=b.setQuery,d.submit=b.submit,d.execute=b.execute,d.$watch("display.query",function(b){d.display.autocomplete=!!b,b&&angular.isString(b)&&a.autocomplete(d.display.query).then(function(a){d.autocomplete.autocompletions=a})}),d.activeIndex=!1,d.selectAutoNext=function(a){a.preventDefault(),d.activeIndex=(d.activeIndex+1)%d.display.queryHistory.length,b.setQuery(d.display.queryHistory[d.activeIndex])},d.selectAutoPrev=function(a){a.preventDefault(),d.activeIndex=(d.activeIndex?d.activeIndex:d.display.queryHistory.length)-1,b.setQuery(d.display.queryHistory[d.activeIndex])},d.fullPageLog=function(){c.path("/terminal"),d.display.hide("log")},d.pathIsTerminal=function(){var a=/^\/terminal.*$/;return a.test(c.path())}}],link:function(b){b.newPage=function(){d.open(d.location.origin,"_blank")},b.newWorkspace=function(){d.open(d.location.origin,"_blank")},b.showMeHow=function(){a.query({name:"Learning Clotho"}).then(function(a){c.path("/trails/"+a[0].id)})},b.goHome=function(){c.path("/")},b.aboutClotho=function(){c.path("/about")},b.teamClotho=function(){c.path("/team")},b.toggleTooltips=function(){console.log("tooltips")}}}}]),angular.module("clotho.commandbar").controller("TerminalCtrl",["$scope","CommandBar",function(a,b){a.log=b.log}]),angular.module("clotho.commandbar").directive("commandBarAutocomplete",["Clotho","CommandBar","$location","$window",function(){return{restrict:"A",templateUrl:"views/_command/autocomplete.html"}}]),angular.module("clotho.commandbar").directive("logEntries",function(){return{restrict:"A",templateUrl:"../views/_command/logEntries.html",scope:{entries:"=logEntries"}}}),angular.module("clotho.commandbar").directive("clothoTokenizer",["$parse",function(){return{restrict:"E",replace:!0,require:"ngModel",scope:{placeholder:"@",startingTags:"=",model:"=ngModel"},templateUrl:"views/_command/tokenizer.html",controller:["$scope","$element","$attrs",function(a){a.tokens=[]}],link:function(a,b,c,d){function e(){console.log("updating model",a.tokens),d.$setViewValue(a.tokens),console.log(d)}a.addToken=function(b){a.tokens.push(b),e()},a.removeToken=function(b){a.tokens.splice(b,1),e()},a.focusInput=function(){b[0].querySelector(".clothoAutocomplete").focus()}}}}]).directive("clothoAutocomplete",["Clotho","$q","$parse","$timeout","$compile","$filter",function(a,b,c,d,e,f){var g=[8,9,13,27,38,40];return{restrict:"A",require:"ngModel",scope:{query:"=ngModel",onSelect:"&autocompleteOnSelect"},controller:["$scope","$element","$attrs",function(){}],link:function(b,c){var h=angular.element("<clotho-autocomplete-listing></clotho-autocomplete-listing>");h.attr({matches:"queryResults",active:"activeIdx",select:"select(activeIdx)",hasFocus:"hasFocus",query:"query"}),b.hasFocus=!1;var i=0,j=function(){b.queryResults=[],b.activeIdx=-1},k=function(c){a.autocomplete(b.query).then(function(a){a&&a.length?b.queryResults=f("limitTo")(a,10):j()})};b.query=void 0;var l;b.$watch("query",function(a){a&&a.length?(b.hasFocus=!0,i>0?(l&&d.cancel(l),l=d(function(){k(a)},i)):k(a)):j()}),b.select=function(a){var e=b.queryResults[a]||b.query;b.onSelect({$item:e}),j(),b.query="",d(function(){c[0].focus()},0,!1)},c.bind("keydown",function(a){-1!==g.indexOf(a.which)&&(a.preventDefault(),8===a.which?b.query.length&&b.$apply(function(){b.query=b.query.substring(0,b.query.length-1)}):40===a.which?(b.activeIdx=(b.activeIdx+1)%b.queryResults.length,b.$digest()):38===a.which?(b.activeIdx=(b.activeIdx?b.activeIdx:b.queryResults.length)-1,b.$digest()):13===a.which||9===a.which?b.$apply(function(){b.select(b.activeIdx)}):27===a.which&&(a.stopPropagation(),j(),b.$digest()))}),c.bind("blur",function(){b.hasFocus=!1}),j(),c.after(e(h)(b))}}}]).directive("clothoAutocompleteListing",function(){return{restrict:"EA",scope:{matches:"=",query:"=",active:"=",hasFocus:"=",select:"&"},replace:!0,templateUrl:"views/_command/autocompleteListing.html",link:function(a){a.isOpen=function(){return a.hasFocus&&a.matches.length>0},a.isActive=function(b){return a.active==b},a.selectActive=function(b){a.active=b},a.selectMatch=function(b){console.log(b),a.select({activeIdx:b})}}}}).directive("clothoAutocompleteMatch",function(){return{restrict:"EA",replace:!0,scope:{index:"=",match:"=",query:"="},templateUrl:"views/_command/autocompleteMatch.html",link:function(){}}}).filter("clothoAutocompleteHighlight",function(){function a(a){return a.replace(/([.?*+^$[\]\\(){}|-])/g,"\\$1")}return function(b,c){return c?b.replace(new RegExp(a(c),"gi"),"<strong>$&</strong>"):b}}).directive("clothoToken",["Clotho",function(a){return{restrict:"E",replace:!0,templateUrl:"views/_command/token.html",require:"ngModel",scope:{model:"=ngModel",onRemove:"&"},controller:["$scope","$element","$attrs",function(){}],link:function(b,c){c.on("click",function(){a.get(b.model.uuid).then(function(){})}),b.removeToken=function(){b.onRemove({model:b.model})}}}}]);
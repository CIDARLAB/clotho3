angular.module("clotho.commandbar",["clotho.core"]),angular.module("clotho.commandbar").service("CommandBar",["Clotho","ClientAPI","$timeout","$q","$document",function(a,b,c,d,e){function f(a){j.unread=j.unread&&!l.log?j.unread+1:1,j.entries.unshift(a),l.show("logSnippet"),j.startLogTimeout()}var g={dateFilter:"mediumTime",timeFilter:"timestamp"},h=function(){return angular.element(e[0].getElementById("clothoCommandBarInput"))},i=function(){return angular.element(e[0].getElementById("clothoCommandBarLogButton"))},j={},k={};k.autocompletions=[],k.autoDetail={},k.detailTemplate={},k.detailModel={},k.detailUUID=-1,j.entries=[{text:"Welcome to Clotho!",from:"server","class":"success",timestamp:Date.now()}];var l={};l.query="",l.queryHistory=[],l.autocomplete=!1,l.autocompleteDetail=!1,l.autocompleteDetailInfo=!1,l.help=!1,l.log=!1,l.logSnippet=!1,l.show=function(a){l[a]||(l[a]=!0)},l.hide=function(a){l[a]&&(l[a]=!1)},l.toggle=function(a){l[a]=!l[a]},l.genLogPos=function(){var a=i()[0];l.logpos={left:a.offsetLeft+a.scrollWidth/2-200+"px",top:a.offsetTop+a.scrollHeight+"px"}},l.genAutocompletePos=function(){var a=h()[0];l.autocompletePos={left:a.offsetLeft+"px",top:a.offsetTop+a.clientHeight+"px"},console.log(a,l)},l.detail=function(b){"undefined"!=typeof b&&(b!=k.detailUUID&&(l.hide("autocompleteDetailInfo"),k.detailUUID=b),a.autocompleteDetail(k.detailUUID).then(function(a){k.autoDetail=a,l.show("autocompleteDetail")}))},l.undetail=function(){l.hide("autocompleteDetail"),l.hide("autocompleteDetailInfo"),k.detailModel={}},l.detailInfo=function(a,b){switch(a){case"command":k.detailTemplate="views/_command/detail-command.html";break;case"author":k.detailTemplate="views/_command/detail-author.html"}k.detailModel=k.autoDetail.sharables[b],"author"==a&&(k.detailModel=k.detailModel.author),l.show("autocompleteDetailInfo")},j.timeout=null,j.startLogTimeout=function(){j.cancelLogTimeout(),j.timeout=c(function(){l.hide("logSnippet")},1e4)},j.cancelLogTimeout=function(){c.cancel(j.timeout)};var m=function(a){console.log("this would be run: "+a),l.hide("autocomplete"),l.undetail()},n=function(c){if(c||(c=l.query),c){var e={"class":"info",from:"client",text:c,timestamp:Date.now()};return l.queryHistory.push(e),l.undetail(),b.say(e),a.submit(c).then(function(a){l.query="",!!a&&b.say({text:a})})}return d.when(!1)};return a.listen("activityLog",function(a){f(a)},"searchbar"),{options:g,display:l,log:j,setQuery:function(a,b){"undefined"!=typeof b&&b.preventDefault(),a&&(l.query=a.value?a.value:a.text)},autocomplete:k,submit:n,execute:m,getCommandBarElement:function(){return elements.commandBarElement},getCommandBarInput:function(){return elements.commandBarInput},focusInput:function(){elements.commandBarInput.focus()}}}]),angular.module("clotho.commandbar").directive("clothoCommandBar",["Clotho","CommandBar","$location","$window",function(a,b,c,d){return{restrict:"A",replace:!0,templateUrl:"views/_command/commandbar.html",controller:["$scope","$element","$attrs",function(d){d.options=b.options,d.log=b.log,d.autocomplete=b.autocomplete,d.display=b.display,d.setQuery=b.setQuery,d.submit=b.submit,d.execute=b.execute,d.$watch("display.query",function(b){d.display.autocomplete=!!b,b&&a.autocomplete(d.display.query).then(function(a){d.autocomplete.autocompletions=a})}),d.prevSubmittedIndex=!1,d.selectAutoNext=function(a){a.preventDefault();var c=d.display.queryHistory.length-1;d.prevSubmittedIndex=d.prevSubmittedIndex?d.prevSubmittedIndex<c?d.prevSubmittedIndex+1:c:0,console.log(d.prevSubmittedIndex),b.setQuery(d.display.queryHistory[d.prevSubmittedIndex])},d.selectAutoPrev=function(a){a.preventDefault();var c=d.display.queryHistory.length-1;d.prevSubmittedIndex=d.prevSubmittedIndex?d.prevSubmittedIndex>0?c:0:d.display.queryHistory.length?c:0,console.log(d.prevSubmittedIndex),b.setQuery(d.display.queryHistory[d.prevSubmittedIndex])},d.fullPageLog=function(){c.path("/terminal"),d.display.hide("log")},d.pathIsTerminal=function(){var a=/^\/terminal.*$/;return a.test(c.path())}}],link:function(b){b.newPage=function(){d.open(d.location.origin,"_blank")},b.newWorkspace=function(){d.open(d.location.origin,"_blank")},b.showMeHow=function(){a.query({name:"Learning Clotho"}).then(function(a){c.path("/trails/"+a[0].id)})},b.goHome=function(){c.path("/")},b.aboutClotho=function(){c.path("/about")},b.teamClotho=function(){c.path("/team")},b.toggleTooltips=function(){console.log("tooltips")}}}}]),angular.module("clotho.commandbar").controller("TerminalCtrl",["$scope","CommandBar",function(a,b){a.log=b.log}]),angular.module("clotho.commandbar").directive("commandBarAutocomplete",["Clotho","CommandBar","$location","$window",function(){return{restrict:"A",templateUrl:"views/_command/autocomplete.html"}}]),angular.module("clotho.commandbar").directive("logEntries",function(){return{restrict:"A",templateUrl:"../views/_command/logEntries.html",scope:{entries:"=logEntries"}}});
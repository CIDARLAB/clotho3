angular.module("clotho.commandbar",["clotho.core"]),angular.module("clotho.commandbar").service("CommandBar",["Clotho","ClientAPI","$timeout","$q","$document",function(a,b,c,d,e){function f(a){l.unread=l.unread&&!n.log?l.unread+1:1,l.entries.unshift(a),console.log(l.entries),n.show("logSnippet"),l.startLogTimeout()}var g={dateFilter:"mediumTime",timeFilter:"timestamp"},h=function(){return angular.element(e[0].getElementById("clothoCommandBar"))},i=function(){return angular.element(e[0].getElementById("clothoCommandBarInput"))},j=function(){return angular.element(e[0].getElementById("clothoCommandBarLogButton"))},k=function(){h().focus()},l={},m={};m.autocompletions=[],m.autoDetail={},m.detailTemplate={},m.detailModel={},m.detailUUID=-1,l.entries=[{text:"Welcome to Clotho!",from:"server","class":"success",timestamp:Date.now()}];var n={};n.query="",n.queryHistory=[],n.autocomplete=!1,n.autocompleteDetail=!1,n.autocompleteDetailInfo=!1,n.help=!1,n.log=!1,n.logSnippet=!1,n.show=function(a){n[a]||(n[a]=!0)},n.hide=function(a){n[a]&&(n[a]=!1)},n.toggle=function(a){n[a]=!n[a]},n.genLogPos=function(){var a=j()[0];n.logpos={left:a.offsetLeft+a.scrollWidth/2-200+"px",top:a.offsetTop+a.scrollHeight+"px"}},n.genAutocompletePos=function(){var a=i()[0];n.autocompletePos={left:a.offsetLeft+"px",top:a.offsetTop+a.clientHeight+"px"}},n.detail=function(b){"undefined"!=typeof b&&(b!=m.detailUUID&&(n.hide("autocompleteDetailInfo"),m.detailUUID=b),a.autocompleteDetail(m.detailUUID).then(function(a){m.autoDetail=a,n.show("autocompleteDetail")}))},n.undetail=function(){n.hide("autocompleteDetail"),n.hide("autocompleteDetailInfo"),m.detailModel={}},n.detailInfo=function(a,b){switch(a){case"command":m.detailTemplate="views/_command/detail-command.html";break;case"author":m.detailTemplate="views/_command/detail-author.html"}m.detailModel=m.autoDetail.sharables[b],"author"==a&&(m.detailModel=m.detailModel.author),n.show("autocompleteDetailInfo")},l.timeout=null,l.startLogTimeout=function(){l.cancelLogTimeout(),l.timeout=c(function(){n.hide("logSnippet")},1e4)},l.cancelLogTimeout=function(){c.cancel(l.timeout)};var o=function(a){console.log("this would be run: "+a),n.hide("autocomplete"),n.undetail()},p=function(c){if(c||(c=n.query),c){var e={"class":"info",from:"client",text:c,timestamp:Date.now()};return n.queryHistory.push(e),n.undetail(),b.say(e),a.submit(c).then(function(a){n.query="",b.say({text:a})})}return d.when(!1)};return a.listen("activityLog",function(a){f(a)},"searchbar"),{options:g,display:n,log:l,setQuery:function(a,b){"undefined"!=typeof b&&b.preventDefault(),a&&(n.query=a.value?a.value:a.text)},autocomplete:m,submit:p,execute:o,getCommandBarElement:h,getCommandBarInput:i,focusInput:k}}]),angular.module("clotho.commandbar").directive("clothoCommandBar",["Clotho","CommandBar","$location","$window",function(a,b,c,d){return{restrict:"A",replace:!0,templateUrl:"views/_command/commandbar.html",controller:["$scope","$element","$attrs",function(d){d.options=b.options,d.log=b.log,d.autocomplete=b.autocomplete,d.display=b.display,d.setQuery=b.setQuery,d.submit=b.submit,d.execute=b.execute,d.$watch("display.query",function(b){d.display.autocomplete=!!b,b&&a.autocomplete(d.display.query).then(function(a){d.autocomplete.autocompletions=a})}),d.prevSubmittedIndex=!1,d.selectAutoNext=function(a){a.preventDefault();var c=d.display.queryHistory.length-1;d.prevSubmittedIndex=d.prevSubmittedIndex?d.prevSubmittedIndex<c?d.prevSubmittedIndex+1:c:0,console.log(d.prevSubmittedIndex),b.setQuery(d.display.queryHistory[d.prevSubmittedIndex])},d.selectAutoPrev=function(a){a.preventDefault();var c=d.display.queryHistory.length-1;d.prevSubmittedIndex=d.prevSubmittedIndex?d.prevSubmittedIndex>0?c:0:d.display.queryHistory.length?c:0,console.log(d.prevSubmittedIndex),b.setQuery(d.display.queryHistory[d.prevSubmittedIndex])},d.fullPageLog=function(){c.path("/terminal"),d.display.hide("log")},d.pathIsTerminal=function(){var a=/^\/terminal.*$/;return a.test(c.path())}}],link:function(b){b.newPage=function(){d.open(d.location.origin,"_blank")},b.newWorkspace=function(){d.open(d.location.origin,"_blank")},b.showMeHow=function(){a.query({name:"Learning Clotho"}).then(function(a){c.path("/trails/"+a[0].id)})},b.goHome=function(){c.path("/")},b.aboutClotho=function(){c.path("/about")},b.teamClotho=function(){c.path("/team")},b.toggleTooltips=function(){console.log("tooltips")}}}}]),angular.module("clotho.commandbar").controller("TerminalCtrl",["$scope","CommandBar",function(a,b){a.log=b.log}]),angular.module("clotho.commandbar").directive("commandBarAutocomplete",["Clotho","CommandBar","$location","$window",function(){return{restrict:"A",templateUrl:"views/_command/autocomplete.html"}}]),angular.module("clotho.commandbar").directive("logEntries",function(){return{restrict:"A",templateUrl:"../views/_command/logEntries.html",scope:{entries:"=logEntries"}}});
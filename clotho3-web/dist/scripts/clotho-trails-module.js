angular.module("clotho.trails",["clotho.foundation","clotho.interface"]),angular.module("clotho.trails").service("Trails",["Clotho","$q",function(a,b){var c=function(c){var d=c.dependencies||null,e=b.defer();if(!d&&!c.mixin)return e.resolve(c),e.promise;var f=[],g=[];return d&&angular.forEach(d,function(b){g.push(a.get(b))}),$clotho.extensions.mixin(c.mixin).then(function(){return b.all(g)}).then(function(a){d={},angular.forEach(a,function(a){d[a.id]=a}),angular.forEach(c.contents,function(a){if("undefined"==typeof a.transclude)f.push(a);else{var b=a.transclude.id,c=a.transclude.chapters;if("all"==c||"undefined"==typeof c)for(var e=0;e<d[b].contents.length;e++)f.push(d[b].contents[e]);else{var g=c.split("-");if(1==g.length)f.push(d[b].contents[g[0]]);else{if(g[0]>g[1])return"wrong format - start must be smaller than end";for(var e=g[0];e<=g[1];e++)f.push(d[b].contents[e])}}}}),c.contents=f,e.resolve(c)}),e.promise},d=function(a,b){var c=b.split("-"),d=a.contents[c[0]].pages[c[1]];return d},e=function(a,b){b="undefined"!=typeof b?b.split("-"):[0,-1];var c;if("undefined"!=typeof a.contents[b[0]].pages[+b[1]+1])c=b[0]+"-"+(+b[1]+1);else{if("undefined"==typeof a.contents[+b[0]+1].pages)return;c=+b[0]+1+"-0"}return c},f=function(a,b){if("0-0"!=b){b="undefined"!=typeof b?b.split("-"):[0,1];var c;return c="undefined"!=typeof a.contents[b[0]].pages[+b[1]-1]?b[0]+"-"+(+b[1]-1):+b[0]-1+"-"+(a.contents[+b[0]-1].pages.length-1)}},g={quiz:"icon-pencil",list:"icon-list-alt",eye:"icon-eye-open",info:"icon-info-sign",video:"icon-film",template:"icon-book",exercise:"icon-edit",undefined:"icon-file"},h=function(a){return a=a||"undefined",g[a]},i=function(a){console.log("favorite trail with id: "+a)};return{compile:c,extractPage:d,calcNextPage:e,calcPrevPage:f,mapIcon:h,share:a.share,favorite:i}}]),angular.module("clotho.trails").directive("trailQuiz",["$http","$templateCache","$compile","Clotho","$interpolate","$q",function(a,b,c,d,e,f){return{restrict:"EA",require:"ngModel",scope:{quiz:"=ngModel",gradeCallback:"=?",advance:"&?"},compile:function(){return{pre:function(d,e){angular.extend(d,d.quiz),d.quiz.question=c("<h5>"+d.quiz.question+"</h5>")(d),a.get("partials/trails/quiz/"+d.quiz.type+"-partial.html",{cache:b}).success(function(a){e.html(c('<div class="quiz">'+a+"</div>")(d))}).error(function(){e.html("<p>Template could not be found...</p>"+JSON.stringify(d.quiz))})},post:function(a){a.createEmptyAnswer=function(b,c){c="undefined"!=typeof c?c:!1,a.quiz.answer=new Array(b.options.length);for(var d=0;d<a.quiz.answer.length;d++)a.quiz.answer[d]=c},a.answerUndefined=function(a){return"undefined"==typeof a.answer||""===a.answer},a.submitQuestion=function(b){d.gradeQuiz(b.questionValue,b.answer,b.answerGenerator).then(function(b){console.log("gradeQuiz result: "+b),a.quiz.submitted=!0,a.quiz.response={},a.quiz.response.result=b,a.gradeCallback(b)})},a.resetQuiz=function(){a.quiz.submitted=!1,a.quiz.response=null,a.quiz.answer=null},a.retryQuiz=function(){if(a.quiz.retry){var b={},c=f.defer();return angular.forEach(a.quiz.retry,function(a,c){b[c]=d.submit(a).then(function(a){return a})}),f.all(b).then(function(b){angular.extend(a.quiz,b),angular.extend(a,a.quiz),a.resetQuiz(),console.log(a),c.resolve()}),c.promise}}}}}}}]),angular.module("clotho.trails").directive("youtube",["Trails","Youtube","$compile","$timeout",function(a,b,c){return{restrict:"EA",replace:!0,scope:{videoId:"@youtube",params:"=?",autoplay:"@?",startMini:"@?",onComplete:"&?"},compile:function(){return{pre:function(){},post:function(a,d){function e(){b.readyPromise.then(function(){f()})}function f(){a.player=new YT.Player(d[0],a.params),$(a.player.a).addClass("youtubePlayer")}var g={width:700,height:525,border:0,autoplay:!1,mini:!1,videoId:a.videoId,playerVars:{autoplay:a.autoplay||a.startMini?1:0,autohide:1,rel:0},events:{}};if(a.params=angular.extend(g,a.params),a.params.videoId){if(a.autoplay=angular.isDefined(a.autoplay)?a.autoplay:a.params.autoplay,a.startMini=angular.isDefined(a.startMini)?a.startMini:a.params.mini,a.params.events.onStateChange=function(b){0==b.data&&a.$apply(a.onComplete())},a.convertToPlayer=function(){f()},a.startMini&&"false"!=a.startMini){a.miniThumb=b.thumbnail(a.videoId,"mqdefault"),b.info(a.videoId).then(function(b){a.miniInfo=b.data,a.miniInfo.durationFormatted=Math.floor(a.miniInfo.duration/60)+":"+(a.miniInfo.duration%60<10?"0":"")+a.miniInfo.duration%60});var h='<div class="row-fluid" style="margin-bottom: 15px"><div class="thumbnail clearfix"><img class="span5" ng-src="{{miniThumb}}"><div class="span7 caption"><h5 style="margin-top:5px">{{ miniInfo.title }}</h5><p style="overflow:hidden; display: -webkit-box; -webkit-line-clamp: 3; -webkit-box-orient: vertical; max-height: 4.5em">{{ miniInfo.description | limitTo:300 }}</p><a class="btn btn-primary" ng-click="convertToPlayer()">Watch Video {{ "(" + miniInfo.durationFormatted +")" }}</a></div></div></div>';d.html(c(h)(a))}else d.html('<img src="../../images/assets/ajax-loader.gif" />'),e();a.$watch("videoId",function(b,c){b!=c&&(a.params=b,f())})}}}}}}]),angular.module("clotho.trails").config(function(){var a=document.createElement("script");a.src="//www.youtube.com/iframe_api";var b=document.getElementsByTagName("script")[0];b.parentNode.insertBefore(a,b)}).service("Youtube",["$http","$rootScope","$q","$timeout",function(a,b,c,d){var e=c.defer();b.$watch(function(){return 1==YT.loaded},function(a){a?e.resolve():d(angular.noop,500)});var f=function(a){var b=/^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;return a.match(b)||a.match(/((\w|-){11})/)?RegExp.$1:!1},g=function(a,b){return b=b||"default","https://img.youtube.com/vi/"+a+"/"+b+".jpg"},h=function(b){return a.get("https://gdata.youtube.com/feeds/api/videos/"+b+"?v=2&prettyprint=true&alt=jsonc").then(function(a){return a.data})};return{readyPromise:e.promise,extract:f,thumbnail:g,info:h}}]),angular.module("clotho.trails").directive("trailPage",["$timeout","$q","$controller",function(a,b,c){return{restrict:"A",template:'<div ng-repeat="comp in pageComponents"><div trail-page-component="comp"></div></div>',scope:{page:"=trailPage"},compile:function(){return{pre:function(d,e){return console.log(d.page),d.setPage=function(a){a&&angular.isObject(a)&&(d.page=a)},d.$watch("page",function(a){a&&d.createPage()}),d.createPage=function(){return d.page.dictionary&&angular.extend(d,d.page.dictionary),$clotho.extensions.css(d.page.css).then(function(){return $clotho.extensions.mixin(d.page.mixin)}).then(function(){return $clotho.extensions.script(d.page.script)}).then(function(){if(d.page.controller){var a={};a.$scope=d;var b=c(d.page.controller,a);e.data("$ngControllerController",b)}}).then(function(){d.pageComponents=d.page.contents},function(a){console.log(a)}).then(function(){return d.page.onload?a(function(){return console.log("loading page onload script"),$clotho.extensions.script(d.page.onload)}):b.when()})},d.createPage(d.page)},post:function(a,b,c){a.$watch(c.trailPage,function(b,c){!!c&&a.setPage(c)},!0)}}}}}]),angular.module("clotho.trails").directive("trailPageComponent",["$compile","$q","$timeout","$http","$templateCache","Youtube",function(a,b,c,d,e,f){return{restrict:"EA",scope:!1,compile:function(){return{pre:function(c,d){var e={};e.hint=function(a){if(!a||!angular.isObject(a))return b.when();var c='<div class="pull-right" hint-button="'+a+'"></div>';return b.when(c)},e.text=function(a){return a?b.when("<div>"+a+"</div>"):b.when()},e.markdown=function(a){return a?b.when("<ui-markdown>"+a+"</ui-markdown>"):b.when()},e.wiki=function(a){return a?b.when("<wiki>"+a+"</wiki>"):b.when()},e.video=function(a){if(!a)return b.when();var d=f.extract(angular.isString(a)?a:a.id);c.videoParams=a.params?a.params:{},c.videoParams.autoplay=angular.isDefined(a.autoplay)?a.autoplay:!1,c.videoParams.mini=angular.isDefined(a.mini)?a.mini:!1;var e='<div><div youtube="'+d+'" params="videoParams"></div></div>';return b.when(e)},e.template=function(a){return a&&angular.isString(a)?(console.log(a),b.when("<div ng-include=\"'"+a+"'\"></div>")):b.when()},e.quiz=function(a){if(!a||!angular.isObject(a))return b.when();if(c.quiz=angular.copy(a),c.quiz.dictionary){var d=angular.copy(c.quiz.dictionary);delete c.quiz.dictionary,angular.extend(c.quiz,d)}c.gradeCallback=function(a){console.log("quiz grade callback result: "+a)};var e='<div trail-quiz ng-model="quiz" grade-callback="gradeCallback()" advance="next()"></div>';return b.when(e)},e.error=function(a){return b.when("<div>"+a+"</div>")};var g=function(a,b){return b&&angular.isString(b)?e[b]?e[b](a):e.error():(console.log("no type passed"),void 0)};c.createPageComponent=function(){return g(c.componentParams,c.componentType).then(function(e){return b.when(d.html(a(e)(c)))})}},post:function(a,b,c){a.$watch(c.trailPageComponent,function(b,c){c&&(a.component=c,a.componentType=a.component.type,a.componentParams=a.component.params,a.createPageComponent())})}}}}}]);
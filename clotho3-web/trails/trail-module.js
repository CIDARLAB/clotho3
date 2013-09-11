'use strict';

//note - requires youtube iFrame API be present
Application.Trails.directive('youtube', ['Trails', '$compile', function(Trails, $compile) {

    return {
        restrict : 'EA',
        replace: true,
        scope: {
            videoId : '@youtube',
            params : '=?',
            autoplay: '@?',
            startMini: '@?',
            onComplete : '&?'
        },
        //template: '<iframe width="{{ params.width }}" height="{{ params.height }}" ng-src="https://www.youtube.com/embed/{{ params.videoId }}?enablejsapi=1&modestbranding=1&rel=0&version=3&playerapiid=trailPlayer&autoplay={{ params.autoplay }}&autohide={{ params.autohide }}&start={{ params.start }}&end={{ params.end }}" frameborder="0" allowfullscreen="1"></iframe>',
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs) {

                },
                post: function postLink(scope, element, attrs) {

                    if (!scope.videoId) return;

                    var videoInfo = Trails.youtubeInfo(scope.videoId);

                    //todo - center if width < 700

                    //defaults
                    var defaults = {
                        width : 700,
                        height : 395,
                        videoId : scope.videoId,
                        playerVars : {
                            autoplay : (!scope.autoplay && !scope.startMini) ? 0 : 1,
                            autohide : 1,
                            rel : 0
                        },
                        events : {}
                    };
                    scope.params = angular.extend(defaults, scope.params);

                    //todo (?) - add more default hooks
                    //todo - avoid overwriting onComplete() if pass in onStateChange()
                    //todo - add $scope.$apply wrapper to all passed in events

                    scope.params.events.onStateChange = function (event) {
                        if (event.data == 0) {
                            scope.$apply(scope.onComplete());
                        }
                    };

                    scope.convertToPlayer = function() {
                        createYoutubePlayer();
                    };

                    if (!!scope.startMini) {
                        scope.miniThumb = Trails.youtubeThumbnail(scope.videoId, 'mqdefault');

                        videoInfo.then(function(json) {
                            scope.miniInfo = json.data;
                            scope.miniInfo.durationFormatted = (Math.floor(scope.miniInfo.duration/60) + ":" + ((scope.miniInfo.duration%60 < 10) ? '0' : '') + (scope.miniInfo.duration%60));
                        });

                        var thumbnailHTML = '<div class="row-fluid thumbnail" style="position:relative;">' +
                            '<img class="span5" ng-src="{{miniThumb}}">' +
                            '<div class="span7 caption">' +
                            '<h5>{{ miniInfo.title }}</h5>' +
                            '<p style="overflow:hidden; display: -webkit-box; -webkit-line-clamp: 3; -webkit-box-orient: vertical; max-height: 4.5em">{{ miniInfo.description | limitTo:300 }}</p>' +
                            '<a class="btn btn-primary" ng-click="convertToPlayer()">Watch Video {{ "(" + miniInfo.durationFormatted +")" }}</a>' +
                            '</div>' +
                            '</div>';

                        element.html($compile(thumbnailHTML)(scope));
                    } else {
                        createYoutubePlayer()
                    }

                    function createYoutubePlayer () {
                        scope.player = new YT.Player(element[0], scope.params);
                    }
                }
            }
        }
    }
}]);


Application.Trails.service('Trails', ['Clotho', '$q', '$dialog', '$http', function(Clotho, $q, $dialog, $http) {

    /**
     * @description Given a URL (youtube.com, youtu.be, watch, embed, etc.), extracts the youtube VideoID. Passing in a VideoId will work.
     * @source http://stackoverflow.com/a/10315969/624466
     * @note assumes length of 11 characters. Youtube may change this in the future.
     *
     * @param {string} url
     * @returns {string} videoId
     */
    var extract_youtube = function(url) {
        var regex = /^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;
        return (url.match(regex) || url.match(/((\w|-){11})/)) ? RegExp.$1 : false;
    };

    //can use youtube names (default, hqdefault, mqdefault, 0, 1, 2, 3)
    //defaults to default (120 x 90)
    var youtubeThumbnail = function(videoId, size) {
        size = size || "default";

        return "https://img.youtube.com/vi/"+ videoId + "/" + size + ".jpg";
    };

    var youtubeInfo = function(videoId) {
        return $http.get('https://gdata.youtube.com/feeds/api/videos/'+videoId+'?v=2&prettyprint=true&alt=jsonc')
            .then(function (data) {
            return data.data
        })
    };


    var compile = function (trail) {

        //If pass by reference (depending on Clotho.get() ) need to copy to don't edit in dependencies
        //trail = angular.copy(trail);

        //todo - better chaining

        var transcludes = trail.dependencies || null,
            deferred = $q.defer();

        if (!transcludes && !trail.mixin) {
            deferred.resolve(trail);
            return deferred.promise;
        }

        var final_contents = [],
            promises = [];

        //get the transcluded trails... will be fast if in collector already
        (transcludes) && angular.forEach(transcludes, function(id) {
            promises.push(Clotho.get(id));
        });

        var mixins = (trail.mixin) ? Application.mixin(trail.mixin) : $q.when();
        mixins
        .then(function() {
            return $q.all(promises)
        })
        //after download all, pluck out the chapters we need
        .then(function (downloads) {

            //reorganize transcludes so can reference by id
            transcludes = {};
            angular.forEach(downloads, function(transclude) {
                transcludes[transclude.id] = transclude;
            });

            //iterate through trail, pushing in chapters
            angular.forEach(trail.contents, function (chapter, ind) {

                if (typeof chapter.transclude == 'undefined') {
                    final_contents.push(chapter);
                } else {
                    //chapters to include :
                    var chapterID = chapter.transclude.id,
                        chapterNum = chapter.transclude.chapters;

                    if ((chapterNum == "all") || (typeof chapterNum == 'undefined')) {
                        for (var i = 0; i < transcludes[chapterID]['contents'].length; i++) {
                            final_contents.push(transcludes[chapterID]['contents'][i]);
                        }
                    } else {
                        var startStop = chapterNum.split("-");
                        if (startStop.length == 1) {
                            final_contents.push(transcludes[chapterID]['contents'][startStop[0]]);
                        } else {
                            if (startStop[0] > startStop[1])
                                return "wrong format - start must be smaller than end";

                            for (var i = startStop[0]; i <= startStop[1]; i++) {
                                final_contents.push(transcludes[chapterID]['contents'][i])
                            }
                        }
                    }
                }
            });

            trail.contents = final_contents;
            deferred.resolve(trail);
        });

        return deferred.promise;
    };

    var favorite = function(id) {
        console.log("favorite trail with id: " + id);
    };

    //todo - better logic...
    var persist = {};

    return {
        extract_youtube : extract_youtube,
        youtubeThumbnail : youtubeThumbnail,
        youtubeInfo : youtubeInfo,
        compile : compile,
        share : Clotho.share,
        favorite : favorite,
        persist : persist
    }
}]);


Application.Trails.controller('TrailMainCtrl', ['$scope', 'Clotho', function($scope, Clotho) {

    $scope.trails = [];

    Clotho.query({schema : "Trail"}).then(function(result) {
        $scope.trails = result;
    });

    $scope.base64icon = base64icon;
}]);

Application.Trails.controller('TrailDetailCtrl', ['$scope', '$route', 'Clotho', 'Trails', '$http', '$timeout', '$templateCache', '$compile', '$keypress', '$q', '$controller', function($scope, $route, Clotho, Trails, $http, $timeout, $templateCache, $compile, $keypress, $q, $controller) {

    //inherited from $routeProvider.resolve clause in application.js
    $scope.id = $route.current.params.id;
    $scope.trail = $route.current.locals.trail;
    $scope.content = $scope.trail.description;

    var load = {};

    load.hint = function loadHint(hint) {
        if (!hint) return $q.when();

        var hintDiv = '<button class="btn pull-right" popover="' + hint + '" popover-trigger="mouseenter" popover-placement="left"><i class="icon-info-sign"></i> Hint</button> ';

        return $q.when(hintDiv);
    };

    load.text = function loadText(text) {
        if (!text) return $q.when();

        //todo - ensure html, if not then surround with tag

        return $q.when(text);
    };

    load.video = function loadVideo(obj) {
        if (!obj) return $q.when();

        var videoId = Trails.extract_youtube( (angular.isString(obj) ? obj : obj.id) );
        $scope.videoParams = (!!obj.params) ? obj.params : {};

        //note - need single outer parent div to compile properly after replace (maybe not in ng-1.2.x)
        //add this attr to move to next automatically on-complete="next()"
        var template = '<div><div youtube="' + videoId + '" params="videoParams" start-mini="'+ (obj.mini || true) +'" autoplay="'+ (obj.autoplay || true) +'"></div></div>';

        //todo - write to avoid timeout? video doesn't update on next() otherwise - probably need to defer instantiation till later
        return $timeout(function() {
           return template;
        });
    };

    load.template = function loadTemplate (url) {
        if (!url) return $q.when();

        //todo - check for url vs. simple template inline

        return $http.get(url, {cache:$templateCache}).then(function (success) {
            return success.data;
        }, function(error) {
            console.log('error retrieving template at: ' + url);
            return '<p class="alert alert-error">That template couldn\'t be found :(</p>'
        });
    };

    load.quiz = function loadQuiz (content) {
        if (!content) return $q.when();

        $scope.quiz = content;
        $scope.gradeCallback = function(data) {
            console.log(data);
        };

        var template = '<div trail-quiz ng-model="quiz" grade-callback="gradeCallback()"></div>';
        return $q.when(template);
    };

    load.error = function loadError(error) {
        return '<h4>Something didn&apos;t work - that type of Page wasn&apos;t recognized</h4>';
    };


    $scope.activate = function(indices) {

        // fields that are handled:
        // backend: CSS (url), mixin (array|url), script (array|url), onload (array|url)
        // content: text (html), video (object), template (url / html), quiz (object)

        //don't activate already active one
        if ($scope.current == indices) return;

        console.log('1 - activating Page');

        $scope.current = indices;
        $scope.content = "";
        //in form <Chapter>-<Page>
        var pos = indices.split("-");
        var page = $scope.trail.contents[pos[0]]['pages'][pos[1]];

        //future in ng-1.2.x, use notify callbacks for updates
        //todo - error callbacks

        console.log(page.intro, page.template, page.video, page.text);

        var loading = $q.defer();
        //$scope.content = '<div class="alert alert-info"><p>Loading content...</p></div>';

        loading.promise = $q.all({
            css : Application.css(page.css),
            mixin: Application.mixin(page.mixin)
        })
        .then(function() {
            //console.log('loading script');

            return Application.script(page.script)

        })
        .then(function (){
            //console.log('loading content');

            if (!!page.dictionary) {
                angular.extend($scope, page.dictionary);
            }

            return $q.all({
                hint : load.hint(page.hint),
                intro : load.text(page.intro),
                video : load.video(page.video),
                template : load.template(page.template),
                quiz : load.quiz(page.quiz),
                outro : load.text(page.outro)
            });
        })
        .then(function(content) {
            console.log(loading, content);

            var contentText = angular.element((content.hint || "") + (content.intro || "") + (content.video || "") + (content.template || "") + (content.quiz || "") + (content.outro || ""));

            //check for controller, must be already included (e.g. by mixin)
            if (page.controller) {
                var locals = {};
                locals.$scope = $scope;
                var ctrl = $controller(page.controller, locals);
                contentText.data('ngControllerController', ctrl);
            }

            $scope.content = $compile(contentText)($scope);

            console.log('loading onload script');
            return Application.script(page.onload);

        }, function(error) {
            //content loading error
            console.log(error);
            //todo - do correctly
            return $q.reject(load.error(error));
        });
    };

    $scope.home = function() {
        $scope.content = $scope.trail.description;
        $scope.current = undefined;
    };

    $scope.favorite = function() {
        //todo - better checking, do initial check
        $scope.favorited = !$scope.favorited;
        Trails.favorite($scope.id);
    };

    $scope.share = function() {
        Trails.share($scope.id)
    };

    $scope.next = function() {
        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, -1],
            newpos;

        if (typeof $scope.trail.contents[oldpos[0]]['pages'][+oldpos[1] + 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
        else if (typeof $scope.trail.contents[+oldpos[0] + 1]['pages'] != 'undefined')
            newpos = (+oldpos[0] + 1) + '-' + 0;
        else {
            $scope.current = undefined;
            return;
        }

        console.log('next is: ' + newpos);

        //todo - force this to run
        $scope.activate(newpos);
    };

    $scope.prev = function() {
        console.log($scope.current);
        if ($scope.current == '0-0') return;

        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, 1],
            newpos;

        if (typeof $scope.trail.contents[oldpos[0]]['pages'][+oldpos[1] - 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
        else if (typeof $scope.trail.contents[+oldpos[0] - 1]['pages'])
            newpos = (+oldpos[0] - 1) + '-' + ($scope.trail.contents[+oldpos[0] - 1]['pages'].length - 1);
        else {
            $scope.current = undefined;
            return;
        }

        //todo - force this to run - as in next()
        $scope.activate(newpos);
    };

    $keypress.on('keydown', {'alt-right' : 'next()', 'alt-left' : 'prev()'}, $scope);

    $scope.base64icon = base64icon;
}]);

Application.Trails.controller('TrailQuizCtrl', ['$scope', 'Clotho', function($scope, Clotho) {

    $scope.createEmptyAnswer = function(quiz, value) {
        value = (typeof value != 'undefined') ? value : false;
        $scope.quiz.answer = new Array(quiz.options.length);
        for (var i = 0; i < $scope.quiz.answer.length; i++) {
            $scope.quiz.answer[i] = value;
        }
    };

    $scope.answerUndefined = function(quiz) {
        return (typeof quiz.answer == 'undefined' || quiz.answer === '');
    };

    $scope.submitQuestion = function(quiz) {
        Clotho.gradeQuiz(quiz).then(function (data) {
            $scope.quiz.submitted = true;
            $scope.quiz.response = data;
        });
    }
}]);

Application.Trails.directive('trailQuiz', ['$http', '$templateCache', '$compile', 'Clotho', function($http, $templateCache, $compile, Clotho) {
    return {
        restrict: "EA",
        require: 'ngModel',
        //could assume function onGrade (or variant) will be present on scope, but I think more maintainable to explicitly declare as attr
        scope: {
            quiz: '=ngModel',
            gradeCallback : '&'
        },

        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs) {

                    $http.get('partials/trails/quiz/' + scope.quiz.type + '-partial.html', {cache: $templateCache})
                        .success(function (data) {
                            element.html($compile(data)(scope));
                        })
                        .error(function(data, status, headers, config) {
                            element.html('<p>Template could not be found...</p>' + JSON.stringify(scope.quiz));
                        });

                },
                post: function postLink(scope, element, attrs) {

                    scope.createEmptyAnswer = function(quiz, value) {
                        value = (typeof value != 'undefined') ? value : false;
                        scope.quiz.answer = new Array(quiz.options.length);
                        for (var i = 0; i < scope.quiz.answer.length; i++) {
                            scope.quiz.answer[i] = value;
                        }
                    };

                    scope.answerUndefined = function(quiz) {
                        return (typeof quiz.answer == 'undefined' || quiz.answer === '');
                    };

                    scope.submitQuestion = function(quiz) {
                        Clotho.gradeQuiz(quiz).then(function (data) {
                            scope.quiz.submitted = true;
                            scope.quiz.response = data;
                            scope.gradeCallback(data);
                        });
                    };

                }
            }
        }
    }
}]);


var base64icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";
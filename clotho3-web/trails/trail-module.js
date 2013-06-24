'use strict';

Application.Trails.service('Trails', ['Clotho', '$q', '$dialog', function(Clotho, $q, $dialog) {

    /**
     * @description Creates a youtube player within an iFrame. Use before createAPIPlayer.
     * @param {object} params as defined:
     * {
     *     videoId:
     *     divId
     *     height:
     *     width:
     *     autoplay:
     *     events: {
     *
     *     }
     *     start:
     *     end:
     *     autoadvance:
    }
    */
    var createFrame = function(params) {
        //if just one parameter, assume its the video id
        if (typeof params == 'string')
            params.videoId = params;

        if (params.cleanId = extract_youtube(params.videoId)) {

            //todo - scale dimensions to max-width:700

            var frame = '<iframe id="'+(params.divId || 'ytplayer')+'" frameborder="0" allowfullscreen="1" ' +
                'width="' + (params.width || 700) + '" ' +
                'height="' + (params.height || 395) + '" ' +
                'src="http://www.youtube.com/embed/'+params.cleanId+'?' +
                'enablejsapi=1&' +
                'modestbranding=1&' +
                'rel=0&' +
                'playerapiid=trailPlayer&' +
                'autoplay=' + (params.autoplay || 1) + '&' +
                'autohide=' + (params.autohide || 1) + '&' +
                (params.start ? 'start=' + params.start : '') + '&' +
                (params.end ? 'end=' + params.end : '') +
                '"></iframe>';

            return frame;
        } else {
            return false;
        }
    };

    var createAPIPlayer = function (params) {
        params = params || {};
        params.events = params.events || {};

        //note: can overwrite certain events if we want to... e.g. params.events.onReady = ...

        /* old code
        var player = new YT.Player(params.divId || 'ytplayer', {
            videoId: cleanId,
            height: params.height || 525,
            width: params.width || 700,
            playervars : {
                enablejsapi: 1,
                modestbranding: 1,
                autohide: params.autohide || 1,
                autoplay: params.autoplay || 1
            },
            events: {
                'onReady': onPlayerReady
                //'onStateChange': onPlayerStateChange
            }
        });
        */

        var player = new YT.Player(document.getElementById(params.divId || 'ytplayer'), {
            events: params.events
        });

        //note: or add custom events...
        /*
        player.addEventListener('onReady', function() {
            console.log("video ready");
        });
        */

        return player;
    };



    /**
     * @description Given a URL (youtube.com, youtu.be, watch, embed, etc.), extracts the youtube VideoID. Passing in a VideoId will work. Adapted from:
     * @url: http://stackoverflow.com/a/10315969/624466
     * @note assumes length of 11 characters. Youtube may change this in the future.
     *
     * @param {string} url
     * @returns {string} videoId
     */
    var extract_youtube = function(url) {
        var regex = /^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;
        return (url.match(regex) || url.match(/((\w|-){11})/)) ? RegExp.$1 : false;
    };


    var compile = function (trail) {

        //need to copy, or compiling will alter the file in the collector (because pass reference)
        trail = angular.copy(trail);

        var transcludes = trail.dependencies || null,
            deferred = $q.defer();

        //if no dependencies listed, don't need to compile
        if (!transcludes) {
            deferred.resolve(trail);
        }

        var final_contents = [],
            promises = [];

        //get the transcluded trails... will be fast if in collector already
        angular.forEach(transcludes, function(uuid) {
            promises.push(Clotho.get(uuid));
        });

        //after download all, pluck out the modules we need
        $q.all(promises).then(function (downloads) {

            //reorganize transcludes so can reference by uuid
            transcludes = {};
            angular.forEach(downloads, function(transclude) {
                transcludes[transclude.uuid] = transclude;
            });

            //iterate through trail, pushing in modules
            angular.forEach(trail.contents, function (mod, ind) {

                //testing
                //console.log(trail.contents[ind]);

                if (typeof mod.transclude == 'undefined') {
                    final_contents.push(mod);
                } else {
                    //modules to include :
                    var modUUID = mod.transclude.uuid,
                        modNum = mod.transclude.modules;

                    if ((modNum == "all") || (typeof modNum == 'undefined')) {
                        for (var i = 0; i < transcludes[modUUID]['contents'].length; i++) {
                            final_contents.push(transcludes[modUUID]['contents'][i]);
                        }
                    } else {
                        var startStop = modNum.split("-");
                        if (startStop.length == 1) {
                            final_contents.push(transcludes[modUUID]['contents'][startStop[0]]);
                        } else {
                            if (startStop[0] > startStop[1])
                                return "wrong format - start must be smaller than end";

                            for (var i = startStop[0]; i <= startStop[1]; i++) {
                                final_contents.push(transcludes[modUUID]['contents'][i])
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

    var share = function(uuid) {
        var dialog_opts = {
            backdrop: true,
            backdropFade: true,
            keyboard: true,
            backdropClick: true,
            templateUrl:  'interface/dialogShare.html',
            controller: 'DialogShareCtrl',
            dependencies : [
                "interface/DialogShareCtrl.js"
            ]
        };

        $dialog.dialog(dialog_opts).open();
    };

    var favorite = function(uuid) {
        console.log("favorite trail with uuid: " + uuid);
    };

    return {
        createFrame : createFrame,
        createAPIPlayer : createAPIPlayer,
        extract_youtube : extract_youtube,
        compile : compile,
        share : share,
        favorite : favorite
    }
}]);


Application.Trails.controller('TrailMainCtrl', ['$scope', 'Clotho', function($scope, Clotho) {
    $scope.trails = Clotho.get('trails');

    $scope.base64icon = base64icon;
}]);

Application.Trails.controller('TrailDetailCtrl', ['$scope', '$route', 'Clotho', 'Trails', '$http', '$timeout', '$templateCache', '$compile', '$keypress', function($scope, $route, Clotho, Trails, $http, $timeout, $templateCache, $compile, $keypress) {

    //inherited from $routeProvider.resolve clause in application.js
    $scope.uuid = $route.current.params.uuid;
    $scope.trail = $route.current.locals.trail;
    $scope.content = $scope.trail.description;

    $scope.loadVideo = function (url) {

        var params = {
            "divId" : "ytplayer",
            "height" : 525,
            "width" : 700,
            "videoId" : url,
            "autoadvance" : true
        };

        $scope.content = Trails.createFrame(params);

        //timeout so content is updated, run after digest()
        $timeout(function() {
            var player = Trails.createAPIPlayer();

            if (params.autoadvance) {
                player.addEventListener('onStateChange', function (event) {
                    if (event.data == 0) {
                        //video is ended, want to advance somehow
                        $scope.next();
                    }
                })
            }
        });
    };

    $scope.loadTemplate = function (url) {
        $http.get(url, {cache:$templateCache}).success(function (data) {
            $scope.content = $compile(data)($scope);
        });
    };

    $scope.loadQuiz = function (content) {
        $scope.quiz = content;

        $http.get('partials/trails/' + content.type + '-partial.html', {cache: $templateCache}).
            success(function (data) {
                $scope.content = $compile(data)($scope);
            });
    };

    $scope.loadPaver = function (paver) {
        var html = '<h4>Something didn&quot;t work - that type of paver wasn&quot;t recognized</h4>';
        $scope.content = html;
    };


    $scope.activate = function(indices) {

        $scope.current = indices;
        $scope.content = "";
        //module, paver
        var pos = indices.split("-");
        var paver = $scope.trail.contents[pos[0]]['pavers'][pos[1]];

        //todo - better handling to ensure script load before next step
        if (paver.script) {
            var paver_name = "paver_" + $scope.uuid + "_" + indices;
            $.getScript(paver.script);
        }

        switch (paver.type) {
            case 'video' : {
                $scope.loadVideo(paver.video);
                break;
            }
            case 'template' : {
                $scope.loadTemplate(paver.template);
                break;
            }
            case 'quiz' : {
                $scope.loadQuiz(paver.content);
                break;
            }
            default : {
                $scope.loadPaver(paver);
            }
        }

    };

    $scope.home = function() {
        $scope.content = $scope.trail.description;
        $scope.current = undefined;
    };

    $scope.favorite = function() {
        //todo - better checking, do initial check
        $scope.favorited = !$scope.favorited;
        Trails.favorite($scope.uuid);
    };

    $scope.share = function() {
        Trails.share($scope.uuid)
    };

    $scope.next = function() {
        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, -1],
            newpos;

        if (typeof $scope.trail.contents[oldpos[0]]['pavers'][+oldpos[1] + 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
        else if (typeof $scope.trail.contents[+oldpos[0] + 1]['pavers'] != 'undefined')
            newpos = (+oldpos[0] + 1) + '-' + 0;
        else {
            $scope.current = undefined;
            return;
        }

        $scope.activate(newpos);
    };

    $scope.prev = function() {
        console.log($scope.current);
        if ($scope.current == '0-0') return;

        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, 1],
            newpos;

        if (typeof $scope.trail.contents[oldpos[0]]['pavers'][+oldpos[1] - 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
        else if (typeof $scope.trail.contents[+oldpos[0] - 1]['pavers'])
            newpos = (+oldpos[0] - 1) + '-' + ($scope.trail.contents[+oldpos[0] - 1]['pavers'].length - 1);
        else {
            $scope.current = undefined;
            return;
        }

        $scope.activate(newpos);
    };

    var x = $keypress.on('keydown', {'right' : 'next()', 'left' : 'prev()'}, $scope);
    console.log(x);

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

var base64icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";
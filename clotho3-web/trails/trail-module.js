'use strict';

Application.Trails.service('Trails', ['Clotho', '$q', function(Clotho, $q) {

    /**
     * @description Creates a youtube player within an iFrame
     * @param {object} params as defined:
     * {
     *     divId
     *     height:
     *     width:
     *     videoId:
     *     events:

     *     //todo
     *     start:
     *     end:
     *     autoplay:
    }
    */
    var createPlayer = function(params) {
        //if just one parameter, assume its the video id
        if (typeof params == 'string')
            params.videoId = params;

        var cleanId = extract_youtube(params.videoId);
        if (cleanId) {
            //todo make sure the div exists, what happens if doesn't? new video if does, not new player
            //see http://stackoverflow.com/questions/7988476/listening-for-youtube-event-in-javascript-or-jquery/7988536#7988536

            new YT.Player(params.divId || 'ytplayer', {
                autohide: 1,
                autoplay: 1,
                height: params.height || 525,
                width: params.width || 700,
                videoId: cleanId,
                events: {
                    //'onReady': onPlayerReady,
                    //'onStateChange': onPlayerStateChange
                }
            });

            return true;
        } else {
            return false;
        }
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


    var compile_trail = function (trail) {

        var transcludes = trail.dependencies || null;

        //if no dependencies listed, don't need to compile
        if (!transcludes) {
            return trail;
        }

        var final_contents = [],
            promises = [];

        //get the transcluded trails
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
        });

        return trail;
    };

    return {
        createPlayer : createPlayer,
        extract_youtube : extract_youtube,
        compile_trail : compile_trail
    }
}]);


Application.Trails.controller('TrailMainCtrl', ['$scope', 'Clotho', function($scope, Clotho) {
    $scope.trails = Clotho.get('trails');

    $scope.base64icon = base64icon;
}]);

Application.Trails.controller('TrailDetailCtrl', ['$scope', '$route', 'Clotho', 'Trails', '$http', function($scope, $route, Clotho, Trails, $http) {

    //todo - move logic into resolve
    Clotho.get($route.current.params.uuid).then(function(result) {
        $scope.trail = Trails.compile_trail(result);
        $scope.content = $scope.trail.description;
    });

    $scope.loadVideo = function (url) {
        //todo - check for start time option

        /*var params = {
            "divId" : "ytplayer",
            "height" : 525,
            "width" : 700,
            "videoId" : url
        };

        Trails.createPlayer(params);*/


        //OLD WAY
        var html = '<iframe id="ytplayer" type="text/html" width="700" height="525" src="'+url+'?autoplay=1&rel=0&autohide=1&enablejsapi=1&version=3&playerapiid=ytplayer" frameborder="0" allowfullscreen></iframe>';
        $scope.content = html;

    };

    $scope.loadTemplate = function (url) {
        //future - once clotho API can return templates, use Clotho.get
        $http.get(url).success(function (data) {
            $scope.content = data;
        });
    };

    $scope.loadQuiz = function (url) {
        $scope.content = '<h4>this would be a quiz once this works</h4>'
    };

    $scope.loadPaver = function (paver) {
        var html = '';
        $scope.content = html;
    };


    $scope.activate = function(indices) {
        //todo - clear video area / update with new video

        //easier to compare
        $scope.current = indices;
        //module, paver
        var pos = indices.split("-");
        var paver = $scope.trail.contents[pos[0]]['pavers'][pos[1]];

        //todo - better handling to ensure script load
        if (paver.script) {
            $script(paver.script);
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
                $scope.loadQuiz(paver.quiz);
                break;
            }
            default : {
                $scope.loadPaver(paver);
            }
        }

    };

    $scope.home = function() {
        $scope.content = $scope.trail.description;
        $scope.current = null;
    };

    $scope.favorite = function() {
        console.log("favorite this trail");
    };

    $scope.share = function() {
        console.log("share this trail");
    };

    $scope.next = function() {
        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, -1],
            newpos;

        if ($scope.trail.contents[oldpos[0]]['pavers'][+oldpos[1] + 1])
            newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
        else if ($scope.trail.contents[+oldpos[0] + 1]['pavers'])
            newpos = (+oldpos[0] + 1) + '-' + 0;
        else {
            $scope.current = null;
            return;
        }

        $scope.activate(newpos);
    };

    $scope.prev = function() {
        var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, -1],
            newpos;

        if ($scope.trail.contents[oldpos[0]]['pavers'][+oldpos[1] - 1])
            newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
        else if ($scope.trail.contents[+oldpos[0] - 1]['pavers'])
            newpos = (+oldpos[0] - 1) + '-' + ($scope.trail.contents[+oldpos[0] - 1]['pavers'].length - 1);
        else {
            $scope.current = null;
            return;
        }

        $scope.activate(newpos);
    };

    $scope.base64icon = base64icon;
}]);

var base64icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";

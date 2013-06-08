'use strict';

Application.Search.directive('clothoSearchbar', ['Clotho', function(Clotho) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl: "search/searchbar-container.html",
        controller: function($scope, $element, $attrs) {
            $scope.query = "";
            $scope.autocompletions = [];
            $scope.autoDetail = {};
            $scope.detailTemplate = {};

            //visibility
            $scope.display = {};
            $scope.display.status = false; //
            $scope.display.autocomplete = false; // autocomplete list
            $scope.display.autocompleteDetail = false; //pane to left of autocomplete
            $scope.display.autocompleteDetailInfo = false; // e.g. command or author
            $scope.display.help = false; // help menu far right
            $scope.display.log = false; // activity log
            $scope.display.logSnippet = false; // snippet right of log button

            //functions

            $scope.$watch('query', function(newValue, oldValue) {
                $scope.display.autocomplete = !!newValue;
                if (!!newValue) {
                    Clotho.autocomplete($scope.query);
                }
            });

            Clotho.listen('autocomplete', function(data) {
                $scope.autocompletions = data;
            }, 'searchbar');

            //pass in $event
            $scope.toggleLog = function(event) {
                if (!$scope.display.logpos) {
                    //get the button
                    var target = event.target;
                    if (target.tagName == "I") {
                        target = target.parentElement;
                    }

                     //calculate where to show
                     $scope.display.logpos = function() {
                         //assumes log is 360px wide
                         return {
                             left : (target.offsetLeft + (target.scrollWidth / 2) - 180) + "px",
                             top : (target.offsetTop + target.scrollHeight)  + "px"
                         };
                     };

                }
                $scope.toggle('log');
            };

            $scope.display.show = function(field) {
                if (!$scope.display[field])
                    $scope.display[field] = true;
            };
            $scope.display.hide = function(field) {
                if ($scope.display[field])
                    $scope.display[field] = false;
            };
            $scope.toggle = function(field) {
                $scope.display[field] = !$scope.display[field];
            };


            $scope.detail = function(uuid) {
                if (typeof uuid == 'undefined') return;

                if (uuid != $scope.detailUUID) {
                    $scope.display.hide('autocompleteDetailInfo');
                    $scope.detailUUID = uuid;
                }

                /*Clotho.get($scope.detailUUID).then(function(result) {
                    $scope.autoDetail = result;
                    $scope.display.show('autocompleteDetail');
                });*/
                Clotho.autocompleteDetail($scope.detailUUID).then(function(result) {
                    $scope.autoDetail = result;
                    $scope.display.show('autocompleteDetail');
                });

            };

            $scope.undetail = function() {
                $scope.display.hide('autocompleteDetail');
                $scope.display.hide('autocompleteDetailInfo');
                $scope.detailModel = {};
            };

            $scope.showAuto = function() {
                if ($scope.query) {
                    $scope.display.show('autocomplete')
                }
            };

            $scope.submit = function() {
                if (!!$scope.query) {
                    Clotho.submit($scope.query);
                    $scope.display.autocomplete = false;
                    $scope.undetail();
                }
            };

            $scope.execute = function(uuid) {
                console.log("this would be run: " + uuid);
                $scope.display.hide('autocomplete');
                $scope.undetail();
            };

            /*** help icons ***/

            $scope.newPage = function() {
                window.open("http://localhost:8000/app/index.html", "_blank");
            };

            $scope.newWorkspace = function() {
                window.open("http://localhost:8000/app/index.html#/trails", "_blank");
            };

            $scope.showMeHow = function() {
                console.log("tutorial");
            };

            $scope.aboutClotho = function() {
                console.log("about clotho");
            };

            $scope.toggleTooltips = function() {
                console.log("tooltips");
            };

            //testing

            $scope.sayTest = function() {
                Clotho.say('This is a test message');
            }


        },
        link: function (scope, element, attrs, controller) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarHelppane', ['Clotho', function(Clotho) {

    return {
        restrict: 'A',
        replace: true,
        template: "",
        controller: function($scope, $element, $attrs) {

        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarAutocomplete', ['Clotho', 'Socket', function(Clotho, Socket) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl: 'search/autocomplete-partial.html',
        controller: function($scope, $element, $attrs) {

            //todo - using index is messy especially if sort
            $scope.detailInfo = function(type, index) {
                //choose template
                switch (type) {
                    case 'command' : {
                        $scope.detailTemplate = 'search/detail-command.html';
                        break;
                    }
                    case 'author' : {
                        $scope.detailTemplate = 'search/detail-author.html';
                        break;
                    }
                    default : {}
                }
                //choose model
                $scope.detailModel = $scope.autoDetail.versions[index];
                if (type == "author")
                    $scope.detailModel = $scope.detailModel.author;

                $scope.display.show('autocompleteDetailInfo');
            };

        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarLog', ['Clotho', '$timeout', function(Clotho, $timeout) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl : "search/log-partial.html",
        controller: function($scope, $element, $attrs) {

            $scope.dateFilter = 'short';
            $scope.timeFilter = 'timestamp';

            $scope.log = {"entries" : [
                {
                    "text" : "Sending message failed",
                    "from" : "client",
                    "class" : "text-error",
                    "timestamp" : 1288399999999
                },
                {
                    "text" : "This is a warning",
                    "from" : "server",
                    "class" : "text-warning",
                    "timestamp" : 1288999999999
                },
                {
                    "text" : "Yay first message worked",
                    "from" : "server",
                    "class" : "text-success",
                    "timestamp" : 1188323623006
                },
                {
                    "text" : "By the way these are automatically sorted by date. This is a really long message to demonstrate what it looks like...",
                    "from" : "client",
                    "class" : "muted",
                    "timestamp" : 1289999908979
                }


            ]};

            $scope.receive = function(data) {
                //todo - get to display for a little bit to the right of the button
                //$scope.display.show('logSnippet');
                //$timeout($scope.display.hide('logSnippet'), 5000);
                $scope.log.entries.unshift(data);
            };

            Clotho.listen("activityLog", function (data) {
                $scope.receive(data);
            }, $scope.$id)

        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);



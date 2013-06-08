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
                $scope.display[field] = true;
            };
            $scope.display.hide = function(field) {
                $scope.display[field] = false;
            };
            $scope.toggle = function(field) {
                $scope.display[field] = !$scope.display[field];
            };


            $scope.detail = function(uuid) {
                Clotho.autocompleteDetail(uuid);
                $scope.detailUUID = uuid;
                $scope.display.show('autocompleteDetail');
                //$scope.autoDetail = uuid;
            };

            $scope.undetail = function(event) {
                $scope.display.hide('autocompleteDetail');
                $scope.display.hide('autocompleteDetailInfo');
                $scope.detailModel = {};
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
                window.open("http://localhost:8080/app/index.html", "_blank");
            };

            $scope.newWorkspace = function() {
                window.open("http://localhost:8080/app/index.html#/trails", "_blank");
            };

            $scope.showMeHow = function() {

            };

            $scope.aboutClotho = function() {

            };

            $scope.toggleTooltips = function() {

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

            /*
             //indexed by uuid - not an array -- could sort in angular if have ranking #
             $scope.autocompletions = {
             "1234567890" : {
             "text" : "This is a command",
             "type" : "command"
             },
             "qwertyuiop" : {
             "text" : "This is a phrase",
             "type" : "phrase"
             },
             "asdfghjk" : {
             "text" : "This is another command",
             "type" : "command"
             },
             "zxcvbnm" : {
             "text" : "Reverse Complement pca1502",
             "type" : "command"
             },
             "xxxxxxxxxxxx" : {
             "text" : "Reverse complement this",
             "type" : "phrase"
             }
             };
             */

            /*$scope.autocompletions = [
                {
                    "uuid" : "1234567890",
                    "text" : "This is a command",
                    "type" : "command"
                },
                {
                    "uuid" : "qwertyuiop",
                    "text" : "This is a phrase",
                    "type" : "phrase"
                },
                {
                    "uuid" : "817924532",
                    "text" : "Reverse Complement pca1502",
                    "type" : "command"
                },
                {
                    "uuid" : "xxxxxxxxxxxx",
                    "text" : "Reverse complement this",
                    "type" : "phrase"
                }
            ];
            */
            $scope.autoDetail = Clotho.get("detail_1234567890");


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


Application.Search.directive('clickOutside', ['Clotho', function(Clotho) {

}]);


//testing
Application.Search.directive('myMouseleave', function() {
    return function(scope, element, attr) {
        element.bind('mouseleave', function(event) {
            //todo - prevent on mouseleave child element
            scope.$eval(attr['myMouseleave'], {$event:event});
        });
    };
});


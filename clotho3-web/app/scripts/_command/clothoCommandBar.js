/**
 * @name clotho-command-bar
 *
 * @description
 * Creates the Clotho Command Bar
 *
 * @example
 * <clotho-command-bar></clotho-command-bar>
 */

angular.module('clotho.commandbar')
.directive('clothoCommandBar', function(Clotho, CommandBar, ClothoCommandHistory, terminalAsideOptions, $location, $window, $compile, $timeout, $clothoModal) {

    return {
        restrict: 'EA',
        scope: true,
        templateUrl: "views/_command/commandbar.html",
        controller: function($scope, $element, $attrs) {

            $scope.logEntries = ClothoCommandHistory.entries;
            $scope.display = CommandBar.display;

            //login

            var showClothoLoginModal = false;
            $scope.toggleLogin = function (force) {
                showClothoLoginModal = angular.isDefined(force) ? force : !showClothoLoginModal;
                if (showClothoLoginModal) {
                    $clothoModal.create({
                        title : 'Clotho Login',
                        'template-url' : "'views/_command/simpleLogin.html'"
                    });
                } else {
                    $clothoModal.destroy();
                }
            };

            //log and snippet

            var logTimeout = null;
            var logTimeoutDelay = 10000;

            $scope.startLogTimeout = function() {
                $scope.cancelLogTimeout();

                logTimeout = $timeout( function() {
                    $scope.display.toggle('logSnippet', false);
                }, logTimeoutDelay);
            };

            $scope.cancelLogTimeout = function() {
                $timeout.cancel(logTimeout);
            };

            //listen for changes to logEntries
            $scope.$watchCollection('logEntries', function () {
                $scope.display.toggle('logSnippet', true);
                $scope.startLogTimeout();
            });

        },
        link : function clothoCommandBarLink(scope, element, attrs, controller) {

            scope.toggleTerminalAside = ClothoCommandHistory.toggleTerminal;

            /*** help icons ***/

            scope.showMeHow = function() {
                Clotho.query({name: 'Learning Clotho'}, {mute: true})
                .then(function (results) {
                    Clotho.startTrail(results[0].id);
                });
            };

            scope.goHome = function() {
                $location.path('/');
            };

            scope.goBrowser = function() {
                $location.path('/browser');
            };

            scope.goEditor = function() {
                $location.path('/editor');
            };

            scope.goTrailBrowser = function() {
                $location.path('/trails');
            }

            scope.aboutClotho = function() {
                $location.path('/about')
            };

            scope.teamClotho = function() {
                $location.path('/team');
            };
        }
    }
});

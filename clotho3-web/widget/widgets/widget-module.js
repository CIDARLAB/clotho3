'use strict';

//todo - incorporate clothoRoot services --- decide how

angular.module('widgetApp', ['clotho.core'])
    //todo - rewrite so still function but not reliant on route
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            otherwise({
                templateUrl:'widget/widgets/widget-partial.html'
            })
    }])
    .run(['$rootScope', 'Clotho', function($rootScope, Clotho) {
        $rootScope.Clotho = Clotho;
    }])
    .directive('widgetDir', function() {
        return {
            restrict: 'A',
            replace: true,
            controller: function($scope, $element) {
                $scope.message = "this is a message from the scope";
            },
            template : '<div class="well">' +
                '<p>This is a directive</p>' +
                'Here is something from the directive controller:' +
                '{{message}}' +
                '</div>'
        }
    })
    .controller('widgetCtrl', ['$scope', 'Clotho', function($scope, Clotho) {
        $scope.greeting = "string from controller";

        Clotho.listen('testMessage', function(data) {
            $scope.greeting = data;
        }, 'widget');
    }]);

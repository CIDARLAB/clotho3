'use strict';

angular.module('widgetApp2', ['clothoPackage'])
    //todo - rewrite so still function but not reliant on route
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            otherwise({
                templateUrl:'widgets/widget2-partial.html'
            })
    }])
    .run(['$rootScope', 'Clotho', function($rootScope, Clotho) {
        $rootScope.Clotho = Clotho;
    }])
    .directive('widgetDir2', function() {
        return {
            restrict: 'A',
            replace: true,
            controller: function($scope, $element) {
                $scope.message = "this is a message from the scope";
            },
            template : '<div class="well">' +
                '<p>This is a directive.... TYPE TWO</p>' +
                'Here is something from the directive controller:' +
                '{{message}}' +
                '</div>'
        }
    })
    .controller('widgetCtrl2', ['$scope', 'Clotho', function($scope, Clotho) {
        $scope.greeting = "string from controller numero dos";

        Clotho.listen('testMessage', function(data) {
            $scope.greeting = data;
        }, 'widget2');
    }]);


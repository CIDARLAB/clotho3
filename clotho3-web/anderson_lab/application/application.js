'use strict';

var Application = Application || {};

Application.Search = angular.module('clotho.search', []);
Application.Interface = angular.module('clotho.interface', []);
Application.Foundation = angular.module('clotho.setup', [])
    .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {

        //extend scope with Clotho API. Don't need to do this in each controller.
        $rootScope.Clotho = Clotho;
    }]);


angular.module('andersonLab', ['clotho.search', 'clotho.setup', 'clotho.interface'])
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider
            .when('/', {
                templateUrl:'partials/home.html'
            })
            //todo - add resolve to query() for people
            .when('/people', {
                templateUrl : 'partials/people.html'
            })
            //todo - add resolve to get() person
            .when('/people/:id', {
                templateUrl:'partials/person.html'
            })
            .when('/research', {
                templateUrl : 'partials/research.html'
            })
            .when('/software', {
                templateUrl : 'partials/software.html'
            })
            .otherwise({
                redirectTo: '/'
            })
    }])
    .run([function () {
        //angular function extensions
        var ext = {};

        ext.isEmpty = function(value) {
            return angular.isUndefined(value) || value === '' || value === null || value !== value;
        };

        ext.isScope = function(obj) {
            return obj && obj.$evalAsync && obj.$watch;
        };

        angular.extend(angular, ext);

    }])
    .controller('PeopleCtrl', ['$scope', function($scope) {

    }])
    .controller('PersonCtrl', ['$scope', '$route', function($scope, $route) {

        //inherit from routeProvider.resolve(), promise fulfilled before render
        $scope.indiv = $route.current.locals.individual;

        $scope.id = $route.current.params.id;

    }])
    .directive('slideshow', [function()]);
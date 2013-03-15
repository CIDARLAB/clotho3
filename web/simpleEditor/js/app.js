'use strict';

/* App Module */

angular.module('simpleEditor', ['editorServices', 'editorDirectives']).
    config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            when('/', {redirectTo:'/first'}).
            when('/:instID', {templateUrl:'partials/interface.html', controller:EditorMainCtrl}).
            when('/:instID/edit', {templateUrl:'partials/interface_edit.html', controller:EditorEditCtrl}).
            otherwise({redirectTo:'/'});
}]);

'use strict';

$clotho.extensions.controller('constructionFiles_wetLabSyntaxCtrl', function ($scope, $http) {

})
//todo - refactor to our popover service?
//todo - CSS
.directive('cfhint', function () {
	return {
		restrict : 'A',
		transclude: true,
		scope : {
			title: '@?cfhintTitle',
			content: '@cfhint'
		},
		template : '<span popover-title="{{title}}" popover-trigger="mouseenter" popover-html-unsafe="{{content}}" popover-append-to-body="true" ng-transclude>'
	}
});
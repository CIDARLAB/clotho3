'use strict';

angular.module('clotho.interface')
.directive('clothoSharable', function ($compile, $http, $templateCache) {
	return {
		restrict: "A",
		scope: {
			content: '=clothoSharable'
		},
		compile: function compile(element, attrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs, controller) {

					//todo - use ClothoSchemas to determine type

					var urlStr = 'views/_interface/sharables/' + 'default' + '.html';

					$http.get(urlStr, {cache: $templateCache})
						.error(function (data, status, headers, config) {
							$http.get('views/_interface/sharables/default.html', {cache: $templateCache})
								.success(function (data, status, headers, config) {
									//force storage for template that didn't exist
									$templateCache.put(urlStr, data);

									element.html(data);
									$compile(element.contents())(scope);
								});
						})
						.success(function (data, status, headers, config) {
							element.html(data);
							$compile(element.contents())(scope);
						});
				},
				post: function () {}
			}
		}
	};
});
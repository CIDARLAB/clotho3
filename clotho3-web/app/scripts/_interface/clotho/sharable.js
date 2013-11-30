'use strict';

angular.module('clotho.interface').directive('sharable', function ($compile, $http, $templateCache) {
	var linker = function (scope, element, attrs) {
		//testing
		scope.base64icon = base64icon;
	};

	return {
		restrict: "E",
		replace: true,
		scope: {
			content: '='
		},
		compile: function compile(element, attrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs, controller) {

					var template;
					if (typeof scope.template != 'undefined')
						template = angular.lowercase(scope.template);
					else {
						template = (scope.content.type != 'Instance') ?
							angular.lowercase(scope.content.type) :
							//testing - temporary hack
							angular.lowercase((/.*\.(.*)$/gi).exec(scope.content.schema.name)[1])
					}

					var urlStr = 'views/_interface/sharables/' + template + '.html';

					$http.get(urlStr, {cache: $templateCache})
						.error(function (data, status, headers, config) {
							//todo - better handling?
							if (status == '404') {
								$http.get('views/_interface/sharables/default.html', {cache: $templateCache})
									.success(function (data, status, headers, config) {
										//force storage for template that didn't exist
										$templateCache.put(urlStr, data);

										element.html(data);
										$compile(element.contents())(scope);
									});
							}
						})
						.success(function (data, status, headers, config) {
							element.html(data);
							$compile(element.contents())(scope);
						});
				},
				post: linker
			}
		}
	};
});
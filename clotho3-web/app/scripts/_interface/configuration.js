angular.module('clotho.interface')
/**
 * simple (temporary) service for handling UI control
 * note that this is not inherited in isolate scopes...you can use $rootScope.interfaceConfig assuming service is exposed on $rootScope
 */
.service('interfaceConfig', function ($injector) {

		return {
			//temporary hack so trails won't show login
			loginVisible : $injector.has('clothoEditorDirective')
		}
});
angular.module('clotho.utils', [])
	.service('ClothoUtils', function() {

		/**
		 * @name verifyUUID
		 * @param uuid
		 * @description Verify that a UUID string is valid
		 */
		var validUUID = function (uuid) {
			return angular.isString(uuid) && uuid.length == 16 && (/[a-zA-Z0-9]{16}/).test(uuid);
		};

		//todo
		var downloadDependencies = function () {};

		return {
			validUUID : validUUID
		}
	});

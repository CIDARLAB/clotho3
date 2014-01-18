// see https://github.com/btford/angular-socket-io/blob/master/socket.spec.js for an io.socket example

describe('clotho.core Socket', function() {

	beforeEach(module('clotho.core'));

	var spy,
		factory;

	beforeEach(inject(function (Socket) {
		spy = jasmine.createSpy('');
		factory = Socket;
	}));

	describe('#emit', function () {

		it('should have an emit function', function() {
			expect(angular.isFunction(factory.emit)).toBe(true);
		})

	});


	describe('#send', function () {

		it('should have a send function', function() {
			expect(angular.isFunction(factory.send)).toBe(true);
		})

	});

});